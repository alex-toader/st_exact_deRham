"""
2D periodic Voronoi complex builder on T^2 = [0,L)^2.

Builds the full DEC complex: d0 (nE x nV), d1 (nF x nE), star0, star1, star2,
edge vectors, Bloch shifts. Uses scipy.spatial.Voronoi with 3x3 replication
for periodicity.

Construction:
  1. Replicate n_cells seeds on 3x3 grid
  2. Compute Voronoi, select ridges by midpoint in [0,L)^2
  3. Identify periodic vertices by wrap + tolerance
  4. Build edges with Bloch shifts
  5. d0: incidence matrix (canonical orientation: lower vertex first)
  6. d1: polygon walk + BFS orientation propagation (d1 d0 = 0)
  7. star0: Delaunay triangle areas (dual 2-cells of vertices)
  8. star1: dual_len / primal_len (Voronoi ridge_points give correct pair)
  9. star2: 1 / face_area (from central copy of Voronoi regions)

Validation (14 tests):
  T1: Euler V-E+F = 0          T8: sum(face_areas) = Area
  T2: nE=3nF, nV=2nF           T9: each edge has 2 faces
  T3: d1 d0 = 0                T10: each face is closed polygon
  T4: star0,star1,star2 > 0    T11: d1 columns are (+1,-1)
  T5: G = Area * I_2           T12: star0 cross-check (independent)
  T6: Betti (1,2,1) on T^2    T13: scalar Laplacian structure (PSD, 1 zero eig)
  T7: sum(star0) = Area        T14: shift consistency (dx = coords[v2]+s*L-coords[v1])

Date: Mar 2026
"""
import numpy as np
from scipy.spatial import Voronoi


def build_2d_periodic_voronoi(n_cells, seed=42, L=1.0, pts=None):
    """
    Build 2D periodic Voronoi DEC complex on T^2 = [0, L)^2.

    Parameters
    ----------
    n_cells : int
        Number of seed points (= number of faces).
    seed : int
        Random seed for reproducibility (ignored if pts is provided).
    L : float
        Box side length.
    pts : ndarray (n_cells, 2), optional
        Explicit seed points in [0, L)^2. If provided, overrides random generation.

    Returns
    -------
    dict with keys:
        d0 : ndarray (nE, nV) — incidence matrix vertices -> edges
        d1 : ndarray (nF, nE) — incidence matrix edges -> faces
        star0 : ndarray (nV,) — vertex dual areas
        star1 : ndarray (nE,) — edge Hodge star = dual_len / primal_len
        star2 : ndarray (nF,) — face Hodge star = 1 / face_area
        dx : ndarray (nE, 2) — edge vectors (in canonical orientation)
        shift : ndarray (nE, 2) int — Bloch shift per edge (units of L)
        coords : ndarray (nV, 2) — canonical vertex coordinates in [0, L)^2
        seeds : ndarray (nF, 2) — seed points (face centers)
        edges : list of (v1, v2) tuples — edge list (canonical orientation)
        faces : list of [v0, v1, ...] lists — ordered vertex lists per face
        nV, nE, nF : int
        L : float
    """
    if pts is None:
        rng = np.random.RandomState(seed)
        pts = rng.rand(n_cells, 2) * L
    else:
        pts = np.asarray(pts, dtype=float)
        assert pts.shape == (n_cells, 2), (
            f"pts shape {pts.shape} != ({n_cells}, 2)")

    # --- Step 1: 3x3 replication for periodicity ---
    rep_pts = []
    rep_cell = []
    central_offset = -1
    for i in range(-1, 2):
        for j in range(-1, 2):
            if i == 0 and j == 0:
                central_offset = len(rep_pts)
            for c in range(n_cells):
                rep_pts.append(pts[c] + np.array([i * L, j * L]))
                rep_cell.append(c)
    rep_pts = np.array(rep_pts)
    rep_cell = np.array(rep_cell)
    assert central_offset >= 0, "central image not found in 3x3 replication"
    vor = Voronoi(rep_pts)

    eps = 1e-10

    # --- Step 2: Select ridges in fundamental domain ---
    # A ridge (= primal edge) is "in" if its midpoint is in [0, L)^2.
    fund_ridges = []
    for ri, (rv1, rv2) in enumerate(vor.ridge_vertices):
        if rv1 == -1 or rv2 == -1:
            continue
        mid = (vor.vertices[rv1] + vor.vertices[rv2]) / 2
        if (-eps <= mid[0] < L - eps) and (-eps <= mid[1] < L - eps):
            fund_ridges.append(ri)

    # --- Step 3: Identify periodic vertices ---
    # Collect all vertex indices used by fundamental ridges.
    used_v = set()
    for ri in fund_ridges:
        rv1, rv2 = vor.ridge_vertices[ri]
        used_v.add(rv1)
        used_v.add(rv2)

    # Group by periodic distance (wrap to [0, L), match within tolerance).
    v2c = {}          # voronoi vertex index -> canonical index
    canon_coords = [] # list of canonical coordinate arrays

    for vi in sorted(used_v):
        v = vor.vertices[vi]
        vw = np.array([v[0] % L, v[1] % L])
        for d in range(2):
            if vw[d] > L - 1e-7:
                vw[d] = 0.0

        found = -1
        for ci, cv in enumerate(canon_coords):
            dv = np.abs(vw - cv)
            dv = np.minimum(dv, L - dv)
            if np.linalg.norm(dv) < 1e-6:
                found = ci
                break

        if found >= 0:
            v2c[vi] = found
        else:
            v2c[vi] = len(canon_coords)
            canon_coords.append(vw.copy())

    nV = len(canon_coords)
    coords = np.array(canon_coords)

    # --- Step 4: Build edges ---
    edge_map = {}  # (c1, c2) with c1 < c2 -> edge index
    edge_data = []

    for ri in fund_ridges:
        rv1, rv2 = vor.ridge_vertices[ri]
        c1, c2 = v2c[rv1], v2c[rv2]
        if c1 == c2:
            continue  # self-loop (degenerate)

        p1, p2 = vor.vertices[rv1], vor.vertices[rv2]
        dx = p2 - p1

        # Bloch shift: periodic image difference
        shift = (np.round((p2 - coords[c2]) / L).astype(int) -
                 np.round((p1 - coords[c1]) / L).astype(int))

        # Canonical orientation: lower vertex index first
        if c1 < c2:
            ekey = (c1, c2)
        else:
            ekey = (c2, c1)
            dx = -dx
            shift = -shift

        if ekey in edge_map:
            continue  # already registered

        # Dual edge (between seed points of the two cells sharing this ridge)
        s1i, s2i = vor.ridge_points[ri]
        dual_len = np.linalg.norm(rep_pts[s2i] - rep_pts[s1i])
        primal_len = np.linalg.norm(dx)
        star1_val = dual_len / primal_len

        edge_map[ekey] = len(edge_data)
        edge_data.append({
            'v1': ekey[0], 'v2': ekey[1],
            'dx': dx, 'shift': shift,
            'star1': star1_val,
            'f1': rep_cell[s1i], 'f2': rep_cell[s2i],
        })

    nE = len(edge_data)
    nF = n_cells

    # --- Step 5: d0 (nE x nV) + edges list ---
    d0 = np.zeros((nE, nV))
    dx_arr = np.zeros((nE, 2))
    shift_arr = np.zeros((nE, 2), dtype=int)
    star1_arr = np.zeros(nE)
    edges = []  # list of (v1, v2) tuples, canonical orientation

    for ei, ed in enumerate(edge_data):
        d0[ei, ed['v1']] = -1
        d0[ei, ed['v2']] = +1
        dx_arr[ei] = ed['dx']
        shift_arr[ei] = ed['shift']
        star1_arr[ei] = ed['star1']
        edges.append((ed['v1'], ed['v2']))

    # --- Step 6: d1 (nF x nE) with consistent orientation ---
    # Build face -> edge adjacency.
    face_edges = {f: [] for f in range(nF)}
    for ei, ed in enumerate(edge_data):
        face_edges[ed['f1']].append(ei)
        face_edges[ed['f2']].append(ei)

    # For each face, walk the polygon boundary and assign consistent signs.
    # Two-step process:
    #   (a) Walk each polygon independently (may have inconsistent global orientation).
    #   (b) BFS propagation to ensure d1[f1,e] + d1[f2,e] = 0 for every edge.
    d1 = np.zeros((nF, nE))
    faces = [[] for _ in range(nF)]  # ordered vertex lists per face

    for f in range(nF):
        fe = face_edges[f]
        if len(fe) == 0:
            continue

        # Build vertex -> edges-of-this-face adjacency.
        v_to_e = {}
        for ei in fe:
            for v in [edge_data[ei]['v1'], edge_data[ei]['v2']]:
                v_to_e.setdefault(v, []).append(ei)

        # Walk: start at first edge, follow vertex chain.
        visited = {fe[0]}
        e0 = fe[0]
        v_start = edge_data[e0]['v1']
        v_curr = edge_data[e0]['v2']
        d1[f, e0] = +1
        faces[f] = [v_start, v_curr]

        for _ in range(len(fe) - 1):
            nxt = [ei for ei in v_to_e[v_curr] if ei not in visited]
            if not nxt:
                break
            en = nxt[0]
            visited.add(en)

            if edge_data[en]['v1'] == v_curr:
                d1[f, en] = +1
                v_curr = edge_data[en]['v2']
            else:
                d1[f, en] = -1
                v_curr = edge_data[en]['v1']
            faces[f].append(v_curr)

        # Remove last vertex if it equals first (closed polygon)
        if len(faces[f]) > 1 and faces[f][-1] == faces[f][0]:
            faces[f] = faces[f][:-1]

        assert len(visited) == len(fe), (
            f"face {f}: walk incomplete ({len(visited)}/{len(fe)} edges)")

    # (b) BFS orientation propagation: ensure each column sums to 0.
    # Face 0 keeps its orientation. For each neighbor face across a shared
    # edge, flip the neighbor's entire row if the shared entry has the same sign.
    from collections import deque
    edge_to_faces = {}
    for e in range(nE):
        edge_to_faces[e] = [f for f in range(nF) if abs(d1[f, e]) > 0.5]

    oriented = np.zeros(nF, dtype=bool)
    oriented[0] = True
    queue = deque([0])

    while queue:
        f = queue.popleft()
        for e in range(nE):
            if abs(d1[f, e]) < 0.5:
                continue
            for g in edge_to_faces[e]:
                if g == f or oriented[g]:
                    continue
                if d1[g, e] * d1[f, e] > 0:
                    d1[g, :] *= -1
                    faces[g] = faces[g][::-1]
                oriented[g] = True
                queue.append(g)

    assert oriented.all(), "BFS did not reach all faces"

    # --- Step 7: star0 (vertex dual areas) ---
    # In 2D Voronoi (z=3): each vertex has 3 incident faces.
    # Dual 2-cell of vertex v = Delaunay triangle connecting the 3 seeds.
    # star0[v] = area of that triangle.
    v_faces = {v: set() for v in range(nV)}
    for ei, ed in enumerate(edge_data):
        v_faces[ed['v1']].add(ed['f1'])
        v_faces[ed['v1']].add(ed['f2'])
        v_faces[ed['v2']].add(ed['f1'])
        v_faces[ed['v2']].add(ed['f2'])

    star0_arr = np.zeros(nV)
    for v in range(nV):
        fl = sorted(v_faces[v])
        if len(fl) != 3:
            raise ValueError(
                "vertex %d has %d incident faces (expected 3 for generic 2D Voronoi). "
                "Mesh has degenerate vertices (4+ cocircular seeds). "
                "Try a different seed." % (v, len(fl)))
        # Get 3 seed positions, closest periodic images to vertex v
        vc = coords[v]
        s = []
        for f in fl:
            sp = pts[f].copy()
            for d in range(2):
                while sp[d] - vc[d] > L / 2:
                    sp[d] -= L
                while sp[d] - vc[d] < -L / 2:
                    sp[d] += L
            s.append(sp)
        # Triangle area (shoelace)
        a, b, c_ = s[0], s[1], s[2]
        star0_arr[v] = 0.5 * abs(
            (b[0] - a[0]) * (c_[1] - a[1]) - (c_[0] - a[0]) * (b[1] - a[1]))

    # --- Step 8: star2 (face inverse areas) ---
    star2_arr = np.zeros(nF)
    face_areas = np.zeros(nF)
    for f in range(nF):
        reg = vor.point_region[central_offset + f]
        region = vor.regions[reg]
        if -1 in region or len(region) == 0:
            raise ValueError(
                f"face {f}: incomplete Voronoi region (vertex at infinity). "
                "3x3 replication should prevent this for central seeds.")
        verts = vor.vertices[region]
        area = 0.5 * abs(sum(
            verts[i, 0] * verts[(i + 1) % len(verts), 1] -
            verts[(i + 1) % len(verts), 0] * verts[i, 1]
            for i in range(len(verts))))
        face_areas[f] = area
        star2_arr[f] = 1.0 / area if area > 0 else 1.0

    return {
        'd0': d0, 'd1': d1,
        'star0': star0_arr, 'star1': star1_arr, 'star2': star2_arr,
        'dx': dx_arr, 'shift': shift_arr,
        'coords': coords, 'seeds': pts,
        'edges': edges, 'faces': faces,
        'nV': nV, 'nE': nE, 'nF': nF, 'L': L,
        'face_areas': face_areas,
    }


# ======================================================================
# Tests
# ======================================================================

def test_complex(n_cells, seed, L=1.0):
    """Run all k=0 validation tests on one complex."""
    c = build_2d_periodic_voronoi(n_cells, seed=seed, L=L)
    d0, d1 = c['d0'], c['d1']
    nV, nE, nF = c['nV'], c['nE'], c['nF']
    s0, s1, s2 = c['star0'], c['star1'], c['star2']
    dx, sh = c['dx'], c['shift']
    area = L ** 2

    errors = []

    # T1: Euler
    euler = nV - nE + nF
    if euler != 0:
        errors.append(f"Euler V-E+F={euler} != 0")

    # T2: Counts
    if nE != 3 * nF:
        errors.append(f"nE={nE} != 3*nF={3*nF}")
    if nV != 2 * nF:
        errors.append(f"nV={nV} != 2*nF={2*nF}")

    # T3: d1 d0 = 0
    dd = d1 @ d0
    dd_err = np.max(np.abs(dd))
    if dd_err > 1e-10:
        errors.append(f"||d1 d0|| = {dd_err:.2e}")

    # T4: Hodge stars positive
    if np.any(s1 <= 0):
        errors.append(f"star1 has non-positive entries: min={s1.min():.4e}")
    if np.any(s2 <= 0):
        errors.append(f"star2 has non-positive entries: min={s2.min():.4e}")
    if np.any(s0 <= 0):
        errors.append(f"star0 has non-positive entries: min={s0.min():.4e}")

    # T5: G = Area * I_2
    G = np.zeros((2, 2))
    for e in range(nE):
        G += s1[e] * np.outer(dx[e], dx[e])
    G_norm = G / area
    diag_dev = max(abs(G_norm[0, 0] - 1), abs(G_norm[1, 1] - 1))
    off_max = max(abs(G_norm[0, 1]), abs(G_norm[1, 0]))
    if diag_dev > 1e-10 or off_max > 1e-10:
        errors.append(f"G/A: diag_dev={diag_dev:.2e}, off={off_max:.2e}")

    # T6: Betti numbers
    rank_d0 = np.linalg.matrix_rank(d0, tol=1e-10)
    rank_d1 = np.linalg.matrix_rank(d1, tol=1e-10)
    b0 = nV - rank_d0
    b1 = (nE - rank_d1) - rank_d0
    b2 = nF - rank_d1
    if (b0, b1, b2) != (1, 2, 1):
        errors.append(f"Betti (b0,b1,b2) = ({b0},{b1},{b2}) != (1,2,1)")

    # T7: sum(star0) = Area (diamond tiling)
    s0_sum = np.sum(s0)
    s0_err = abs(s0_sum / area - 1)
    if s0_err > 1e-8:
        errors.append(f"sum(star0)/Area = {s0_sum/area:.10f} != 1")

    # T8: sum(face_areas) = Area
    fa_sum = np.sum(c['face_areas'])
    fa_err = abs(fa_sum / area - 1)
    if fa_err > 1e-8:
        errors.append(f"sum(face_areas)/Area = {fa_sum/area:.10f} != 1")

    # T9: Each edge has exactly 2 faces
    for e in range(nE):
        n_faces = np.count_nonzero(d1[:, e])
        if n_faces != 2:
            errors.append(f"edge {e} has {n_faces} faces (expected 2)")
            break

    # T10: Each face is a closed polygon (>= 3 edges, forms a cycle)
    for f in range(nF):
        f_edges = np.where(np.abs(d1[f]) > 0.5)[0]
        if len(f_edges) < 3:
            errors.append(f"face {f} has {len(f_edges)} edges (< 3)")
            break
        # Verify cycle: walk the polygon and check we return to start
        v_to_e = {}
        for ei in f_edges:
            for v in range(nV):
                if abs(d0[ei, v]) > 0.5:
                    v_to_e.setdefault(v, []).append(ei)
        # Each vertex in this face should have exactly 2 edges
        bad_v = [v for v, el in v_to_e.items() if len(el) != 2]
        if bad_v:
            errors.append(f"face {f}: vertices {bad_v[:3]} don't have 2 edges")
            break

    # T11: d1 orientation — each column has exactly one +1 and one -1
    for e in range(nE):
        col = d1[:, e]
        nz = col[np.abs(col) > 0.5]
        if len(nz) != 2 or not (
                (nz[0] > 0 and nz[1] < 0) or (nz[0] < 0 and nz[1] > 0)):
            errors.append(f"edge {e}: d1 column not (+1,-1), got {nz}")
            break

    # T12: star0 cross-check — recompute from vertex->face adjacency
    v_faces_check = {v: set() for v in range(nV)}
    for e in range(nE):
        v1e = np.where(d0[e] < -0.5)[0]
        v2e = np.where(d0[e] > 0.5)[0]
        if len(v1e) == 0 or len(v2e) == 0:
            continue
        faces_e = np.where(np.abs(d1[:, e]) > 0.5)[0]
        for ff in faces_e:
            v_faces_check[v1e[0]].add(ff)
            v_faces_check[v2e[0]].add(ff)
    star0_check = np.zeros(nV)
    seeds = c['seeds']
    vcoords = c['coords']
    for v in range(nV):
        fl = sorted(v_faces_check[v])
        if len(fl) != 3:
            errors.append(f"star0 check: vertex {v} has {len(fl)} faces")
            break
        vc = vcoords[v]
        sv = []
        for ff in fl:
            sp = seeds[ff].copy()
            for dd in range(2):
                while sp[dd] - vc[dd] > L / 2:
                    sp[dd] -= L
                while sp[dd] - vc[dd] < -L / 2:
                    sp[dd] += L
            sv.append(sp)
        a, b, cc = sv[0], sv[1], sv[2]
        star0_check[v] = 0.5 * abs(
            (b[0] - a[0]) * (cc[1] - a[1]) - (cc[0] - a[0]) * (b[1] - a[1]))
    s0_match = np.max(np.abs(s0 - star0_check))
    if s0_match > 1e-12:
        errors.append(f"star0 cross-check: max deviation = {s0_match:.2e}")

    # T13: scalar Laplacian structure check
    # Δ₀ = S0^{-1/2} d0^T S1 d0 S0^{-1/2} (symmetric form).
    # Must be: PSD, exactly 1 zero eigenvalue (b0=1), λ1 in correct ballpark.
    S0isq = np.diag(1.0 / np.sqrt(s0))
    S1mat = np.diag(s1)
    lap0 = S0isq @ d0.T @ S1mat @ d0 @ S0isq
    eigs_lap = np.sort(np.linalg.eigvalsh(lap0))
    n_neg = np.sum(eigs_lap < -1e-10)
    n_zero = np.sum(np.abs(eigs_lap) < 1e-8)
    if n_neg > 0:
        errors.append(f"scalar Laplacian: {n_neg} negative eigenvalues (not PSD)")
    if n_zero != 1:
        errors.append(f"scalar Laplacian: {n_zero} zero eigenvalues (expected 1)")
    eig_nz = eigs_lap[np.abs(eigs_lap) > 1e-8]
    expected = (2 * np.pi / L) ** 2
    if len(eig_nz) > 0:
        lap_err = abs(eig_nz[0] / expected - 1)
    else:
        lap_err = 1.0
    # Random Voronoi: large discretization error on coarse meshes.
    # Structural test: λ1 within 50% of continuum value.
    if lap_err > 0.50:
        errors.append(f"scalar Laplacian: λ1/(2π/L)² = {eig_nz[0]/expected:.4f}, "
                       f"err={lap_err:.2e} > 0.50")

    # T14: Bloch shift consistency — dx = coords[v2] + shift*L - coords[v1]
    vcoords2 = c['coords']
    for e in range(nE):
        v1e = np.where(d0[e] < -0.5)[0][0]
        v2e = np.where(d0[e] > 0.5)[0][0]
        dx_expected = vcoords2[v2e] + sh[e] * L - vcoords2[v1e]
        if np.linalg.norm(dx[e] - dx_expected) > 1e-10:
            errors.append(f"edge {e}: shift inconsistent with dx")
            break

    # T15: faces consistent with d1 orientation
    edge_set_bi = set(c['edges']) | {(j, i) for i, j in c['edges']}
    for f in range(nF):
        face = c['faces'][f]
        n_edges_f = np.count_nonzero(np.abs(d1[f]) > 0.5)
        if len(face) != n_edges_f:
            errors.append(f"face {f}: len(faces)={len(face)} != n_edges={n_edges_f}")
            break
        for v_pos in range(len(face)):
            i = face[v_pos]
            j = face[(v_pos + 1) % len(face)]
            if (i, j) not in edge_set_bi:
                errors.append(f"face {f}: ({i},{j}) not an edge")
                break
        else:
            continue
        break

    # T16: faces orientation matches d1
    edge_map_t = {}
    for idx, (i, j) in enumerate(c['edges']):
        edge_map_t[(i, j)] = (idx, +1)
        edge_map_t[(j, i)] = (idx, -1)
    for f in range(nF):
        face = c['faces'][f]
        for v_pos in range(len(face)):
            i = face[v_pos]
            j = face[(v_pos + 1) % len(face)]
            e_idx, orient = edge_map_t[(i, j)]
            if abs(d1[f, e_idx] - orient) > 0.5:
                errors.append(
                    f"face {f}: orientation mismatch at edge {e_idx} "
                    f"(d1={d1[f,e_idx]:.0f}, walk={orient})")
                break
        else:
            continue
        break

    # Report
    status = "PASS" if not errors else "FAIL"
    lap_ratio = eig_nz[0] / expected if len(eig_nz) > 0 else 0
    print(f"  n={n_cells:3d} seed={seed:4d} | nV={nV} nE={nE} nF={nF} | "
          f"V-E+F={euler} | d1d0={dd_err:.0e} | "
          f"G/A diag_dev={diag_dev:.0e} off={off_max:.0e} | "
          f"b=({b0},{b1},{b2}) | s0/A={s0_sum/area:.6f} | "
          f"lam1={lap_ratio:.4f} | {status}")
    for err in errors:
        print(f"    ERROR: {err}")

    return len(errors) == 0, c


if __name__ == '__main__':
    print("=== 2D periodic Voronoi builder validation ===\n")

    configs = [
        (30, 42), (80, 42), (30, 137),
        (30, 999), (50, 2024), (100, 7),
    ]

    all_pass = True
    for n_cells, seed in configs:
        ok, _ = test_complex(n_cells, seed)
        all_pass = all_pass and ok

    print(f"\n{'ALL PASS' if all_pass else 'SOME FAILED'}")
