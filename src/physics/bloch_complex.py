"""
Bloch-periodic cochain complex: d₂(k) construction.

Extends the exactness-preserving recurrence from level 1→2 (d₁, in gauge_bloch.py)
to level 2→3 (d₂). The same pattern applies: derive Bloch phases from the lower
operator (d₁(k)) via BFS on the boundary of each cell.

Recurrence formula:
    d₂[c, f_next] = -d₂[c, f_curr] · d₁[f_curr, e] / d₁[f_next, e]
where e is the shared edge between f_curr and f_next on cell c's boundary.
"""

import numpy as np


def build_face_edge_map(faces, edges):
    """Map face index → set of edge indices (unoriented).

    Args:
        faces: list of face vertex lists
        edges: list of (v1, v2) edge tuples

    Returns:
        dict[int, set[int]]: face_idx → set of edge indices
    """
    E_canonical = [(min(a, b), max(a, b)) for a, b in edges]
    edge_to_idx = {e: i for i, e in enumerate(E_canonical)}
    result = {}
    for f_idx, face in enumerate(faces):
        n = len(face)
        fe = set()
        for i in range(n):
            a, b = face[i], face[(i + 1) % n]
            key = (min(a, b), max(a, b))
            if key in edge_to_idx:
                fe.add(edge_to_idx[key])
        result[f_idx] = fe
    return result


def build_d2_top_from_foam(data):
    """Build topological d₂ from foam data (face_to_cells).

    Determines face orientations relative to cells. Validates at k=0
    that d₂·d₁ = 0.

    Args:
        data: foam dict with V, E, F, L_vec, face_to_cells, cell_centers

    Returns:
        (d2_top, cell_face_incidence):
          d2_top: (nC, nF) real matrix
          cell_face_incidence: list of [(face_idx, ±1), ...] per cell
    """
    from physics.gauge_bloch import compute_edge_shifts, build_d0_bloch, build_d1_bloch_exact

    V, E, F = data['V'], data['E'], data['F']
    L_vec = data['L_vec']
    face_to_cells = data['face_to_cells']
    cell_centers = data['cell_centers']
    nF = len(F)
    nC = len(cell_centers)

    # Initial guess: d₂[cell_a, f] = +1, d₂[cell_b, f] = -1
    d2_top = np.zeros((nC, nF))
    for f_idx in range(nF):
        ca, cb = face_to_cells[f_idx]
        d2_top[ca, f_idx] = +1
        d2_top[cb, f_idx] = -1

    # Validate at k=0
    shifts = compute_edge_shifts(V, E, L_vec)
    k0 = np.zeros(3)
    d0_k0 = build_d0_bloch(V, E, k0, L_vec, shifts)
    d1_k0 = build_d1_bloch_exact(V, E, F, k0, L_vec, d0_k0)

    if np.linalg.norm(d2_top @ d1_k0) > 1e-10:
        # Fix orientations using face normals
        for f_idx in range(nF):
            face = F[f_idx]
            verts = np.array([V[v] for v in face[:3]])
            v1 = verts[1] - verts[0]
            v2 = verts[2] - verts[0]
            v1 -= np.round(v1 / L_vec) * L_vec
            v2 -= np.round(v2 / L_vec) * L_vec
            normal = np.cross(v1, v2)
            if np.linalg.norm(normal) < 1e-10:
                continue
            normal = normal / np.linalg.norm(normal)

            ca, cb = face_to_cells[f_idx]
            centroid = np.mean(np.array([V[v] for v in face]), axis=0)
            dir_a = cell_centers[ca] - centroid
            # Minimum image convention for periodic boundaries
            dir_a -= np.round(dir_a / L_vec) * L_vec

            if np.dot(normal, dir_a) > 0:
                d2_top[ca, f_idx] = -1
                d2_top[cb, f_idx] = +1
            else:
                d2_top[ca, f_idx] = +1
                d2_top[cb, f_idx] = -1

    norm_check = np.linalg.norm(d2_top @ d1_k0)
    assert norm_check < 1e-10, \
        f"d2_top @ d1(k=0) != 0 after orientation fix: {norm_check:.2e}"

    # Build cell_face_incidence list
    cfi = [[] for _ in range(nC)]
    for c in range(nC):
        for f in range(nF):
            if abs(d2_top[c, f]) > 0.5:
                cfi[c].append((f, int(d2_top[c, f])))

    return d2_top, cfi


def build_d2_bloch_exact(cell_face_inc, face_edge_map, d1_k, n_cells, n_faces):
    """Build exactness-preserving d₂(k) via recurrence on cell boundaries.

    For each cell c, BFS over its boundary faces. At each edge shared by
    two faces of c, the constraint d₂d₁=0 determines the phase ratio:
        d₂[c, f_next] = -d₂[c, f_curr] · d₁[f_curr, e] / d₁[f_next, e]

    This is the level 2→3 analogue of the d₁ recurrence at level 1→2.

    Args:
        cell_face_inc: list of [(face_idx, orientation), ...] per cell
        face_edge_map: dict face_idx → set of edge indices
        d1_k: (nF, nE) complex matrix (exact d₁ at wave vector k)
        n_cells: number of cells
        n_faces: number of faces

    Returns:
        d2: (n_cells, n_faces) complex matrix with d₂(k)d₁(k) = 0
    """
    d2 = np.zeros((n_cells, n_faces), dtype=complex)

    for c_idx in range(n_cells):
        faces_of_c = cell_face_inc[c_idx]
        if not faces_of_c:
            continue

        # Build face adjacency within this cell (faces sharing an edge)
        adj = {f: [] for f, _ in faces_of_c}
        for i, (f1, _) in enumerate(faces_of_c):
            for j, (f2, _) in enumerate(faces_of_c):
                if i >= j:
                    continue
                shared = face_edge_map[f1] & face_edge_map[f2]
                for e in shared:
                    adj[f1].append((f2, e))
                    adj[f2].append((f1, e))

        # BFS recurrence from seed face
        f0, o0 = faces_of_c[0]
        d2[c_idx, f0] = o0  # seed: topological orientation
        visited = {f0}
        queue = [f0]
        bfs_edges = set()

        while queue:
            fc = queue.pop(0)
            for fn, e_shared in adj[fc]:
                if fn in visited:
                    continue
                # d₁[f, e] includes face-edge orientation sign.
                # Constraint d₂d₁=0 at edge e: d₂[c,fc]·d₁[fc,e] + d₂[c,fn]·d₁[fn,e] = 0
                d1c = d1_k[fc, e_shared]
                d1n = d1_k[fn, e_shared]
                if abs(d1n) < 1e-14:
                    raise ValueError(
                        f"d1[{fn},{e_shared}] = 0, cannot propagate "
                        f"(cell {c_idx})")
                d2[c_idx, fn] = -d2[c_idx, fc] * d1c / d1n
                visited.add(fn)
                queue.append(fn)
                bfs_edges.add((fc, fn, e_shared))
                bfs_edges.add((fn, fc, e_shared))

        assert len(visited) == len(faces_of_c), \
            f"Cell {c_idx}: BFS visited {len(visited)}/{len(faces_of_c)} faces"

        # Cell holonomy check: non-BFS edges must be consistent
        for fc in adj:
            for fn, e_shared in adj[fc]:
                if (fc, fn, e_shared) in bfs_edges:
                    continue
                if fc >= fn:
                    continue
                d1c = d1_k[fc, e_shared]
                d1n = d1_k[fn, e_shared]
                if abs(d1n) < 1e-14:
                    continue
                expected = -d2[c_idx, fc] * d1c / d1n
                actual = d2[c_idx, fn]
                assert abs(expected - actual) < 1e-10, \
                    f"Cell {c_idx}: holonomy failure on non-BFS edge " \
                    f"({fc},{fn},e={e_shared}): " \
                    f"expected {expected:.6f}, got {actual:.6f}"

    return d2
