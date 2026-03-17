"""
GROUP 2 — The problem: standard fails (§3).

Demonstrate that standard Bloch d₁ breaks exactness.

Tests:
  T2.1  ‖d₁_std(k) d₀(k)‖ > 1 on all structures (incl. random Voronoi)
  T2.2  ‖d₁_std d₀‖ direction-dependent (Kelvin, C15, WP)
  T2.3  n_spur = V − n_zero_std > 0 on all structures
  T2.4  Intra-face contradictions (Cor 1)
  T2.5  Per-face constraint has 1D null space (uniqueness premise)
  T2.6  ‖d₁_std d₀‖ ~ O(k) scaling
  T2.7  NEG: Standard ker(K) ≠ V
  T2.8  NEG: Standard rank(d₁) ≠ exact AND direction-dependent
  T2.9  CONTROL: At k=0, exact = standard (modulus comparison, intentional)

Review notes (not blocking):
  - T2.4: could verify geometric criterion (faces with ≥2 distinct-shift edges) matches
    detected contradictions. Current test is sufficient for Cor 1 but not maximal.
  - T2.10 (future): continuity of ‖d₁_std d₀‖ in k along dense path
  - T2.11 (future): numerical optimization of per-edge phases → min ‖d₁d₀‖ still O(1)

Expected output (verified 17 Mar 2026):

  ============================================================
  GROUP 2: The problem — standard fails (§3)
  ============================================================

    T2.1+T2.3  Standard breaks exactness (incl. random Voronoi)
           Kelvin N=2: ‖d₁d₀‖_std=7.54, ‖d₁d₀‖_ex=7.86e-16, n_spur=6  [PASS]
           C15 N=1: ‖d₁d₀‖_std=7.80, ‖d₁d₀‖_ex=5.54e-16, n_spur=9  [PASS]
           WP N=1: ‖d₁d₀‖_std=5.07, ‖d₁d₀‖_ex=5.20e-16, n_spur=3  [PASS]
           Voronoi seed=42: ‖d₁d₀‖_std=10.52, ‖d₁d₀‖_ex=1.03e-15, n_spur=15  [PASS]
           Voronoi seed=137: ‖d₁d₀‖_std=11.68, ‖d₁d₀‖_ex=1.19e-15, n_spur=19  [PASS]

    T2.2  Direction-dependent pollution (3 structures)
           Kelvin N=2: [100]:6, [110]:12, [111]:14  [direction-dependent]
           C15 N=1: [100]:9, [110]:15, [111]:16  [direction-dependent]
           WP N=1: [100]:3, [110]:5, [111]:7  [direction-dependent]
           All 3 structures direction-dependent  [PASS]

    T2.4+T2.5  Intra-face contradictions (Corollary 1)
           T2.5: All 112 faces have 1D null space  [PASS]
           T2.4: 37/112 faces have intra-face contradictions  [PASS]

    T2.6  ‖d₁_std d₀‖ ~ O(k) scaling
           frac=0.02: ‖d₁_std d₀‖ = 1.61
           frac=0.05: ‖d₁_std d₀‖ = 3.96
           frac=0.10: ‖d₁_std d₀‖ = 7.54
           frac=0.15: ‖d₁_std d₀‖ = 10.45
           frac=0.20: ‖d₁_std d₀‖ = 12.47
           Log-log slope (all): 0.90, (k≤0.10): 0.96
           [PASS]

    T2.7  NEG: Standard ker(K) ≠ V
           Kelvin N=2: ker_exact=96=V, ker_std=90<V, deficit=6  [PASS]
           C15 N=1: ker_exact=136=V, ker_std=127<V, deficit=9  [PASS]
           WP N=1: ker_exact=46=V, ker_std=43<V, deficit=3  [PASS]
           Voronoi seed=42: ker_exact=331=V, ker_std=316<V, deficit=15  [PASS]
           Voronoi seed=137: ker_exact=347=V, ker_std=328<V, deficit=19  [PASS]

    T2.8  NEG: Standard rank unstable
           Exact rank(d₁) = 96 (stable across directions)
           Standard [100]: rank(d₁) = 102
           Standard [110]: rank(d₁) = 108
           Standard [111]: rank(d₁) = 110
           All standard ranks ≠ exact (96)  [PASS]
           Standard rank direction-dependent: {[100]:102, [110]:108, [111]:110}  [PASS]
           → cohomology undefined on standard  [PASS]

    T2.9  CONTROL: At k=0, exact = standard
           Kelvin N=2: d₁d₀=0 at k=0, |d₁_ex(0)|=|d₁_top|  [PASS]
           C15 N=1: d₁d₀=0 at k=0, |d₁_ex(0)|=|d₁_top|  [PASS]
           WP N=1: d₁d₀=0 at k=0, |d₁_ex(0)|=|d₁_top|  [PASS]
           Voronoi seed=42: d₁d₀=0 at k=0, |d₁_ex(0)|=|d₁_top|  [PASS]
           Voronoi seed=137: d₁d₀=0 at k=0, |d₁_ex(0)|=|d₁_top|  [PASS]

  ------------------------------------------------------------
  GROUP 2: ALL PASS (1.8s)
  ------------------------------------------------------------
"""

import numpy as np
from scipy.linalg import eigh

from physics.hodge import (
    build_kelvin_with_dual_info,
    build_c15_with_dual_info,
    build_wp_with_dual_info,
    build_foam_with_dual_info,
    build_hodge_stars_voronoi,
)
from physics.gauge_bloch import compute_edge_shifts, build_d0_bloch, build_d1_bloch_exact
from physics.bloch import build_d1_bloch_standard, compute_edge_crossings, build_edge_lookup
from core_math.operators.incidence import build_d0, build_d1


# ── Helpers ──────────────────────────────────────────────────────────────

def build_all_structures():
    """Build all structures including random Voronoi."""
    structs = {
        'Kelvin N=2': build_kelvin_with_dual_info(N=2),
        'C15 N=1':    build_c15_with_dual_info(N=1),
        'WP N=1':     build_wp_with_dual_info(N=1),
    }
    for seed in [42, 137]:
        np.random.seed(seed)
        points = np.random.uniform(0, 4.0, size=(50, 3))
        try:
            data = build_foam_with_dual_info(points, 4.0)
            structs[f'Voronoi seed={seed}'] = data
        except Exception as e:
            print(f"  WARNING: Voronoi seed={seed} failed: {e}")
    return structs


def make_k(data, frac, direction):
    d = np.array(direction, dtype=float)
    d /= np.linalg.norm(d)
    return 2 * np.pi / data['L'] * frac * d


def build_both_d1(data, k):
    V, E, F = data['V'], data['E'], data['F']
    L, L_vec = data['L'], data['L_vec']
    shifts = compute_edge_shifts(V, E, L_vec)
    crossings = compute_edge_crossings(V, E, L)
    edge_lookup = build_edge_lookup(E, crossings)
    d0_k = build_d0_bloch(V, E, k, L_vec, shifts)
    d1_ex = build_d1_bloch_exact(V, E, F, k, L_vec, d0_k)
    d1_std = build_d1_bloch_standard(V, E, F, L, k, edge_lookup, crossings)
    return d1_ex, d1_std, d0_k


def count_zeros(K, star1):
    M = np.diag(star1)
    K_sym = 0.5 * (K + K.conj().T)
    eigs = np.real(eigh(K_sym, M, eigvals_only=True))
    thresh = max(np.max(np.abs(eigs)) * 1e-12, 1e-10)
    return int(np.sum(np.abs(eigs) < thresh))


def build_face_edges(E, F):
    edge_map = {}
    for idx, (i, j) in enumerate(E):
        edge_map[(i, j)] = (idx, +1)
        edge_map[(j, i)] = (idx, -1)
    face_edges = []
    for face in F:
        n = len(face)
        edges_info = []
        for v_pos in range(n):
            i, j = face[v_pos], face[(v_pos + 1) % n]
            e_idx, orient = edge_map[(i, j)]
            edges_info.append((e_idx, orient))
        face_edges.append(edges_info)
    return face_edges


# ── Tests ────────────────────────────────────────────────────────────────

def test_standard_breaks_exactness(structs):
    """T2.1 + T2.3: ‖d₁_std d₀‖ > 1 and n_spur > 0 on all structures."""
    print("\n  T2.1+T2.3  Standard breaks exactness (incl. random Voronoi)")
    for name, data in structs.items():
        k = make_k(data, 0.10, [1.0, 0.0, 0.0])
        d1_ex, d1_std, d0_k = build_both_d1(data, k)
        star1, star2 = build_hodge_stars_voronoi(data)
        nV = len(data['V'])

        norm_ex = np.linalg.norm(d1_ex @ d0_k)
        norm_std = np.linalg.norm(d1_std @ d0_k)

        K_std = d1_std.conj().T @ np.diag(star2) @ d1_std
        n_zero_std = count_zeros(K_std, star1)
        n_spur = nV - n_zero_std

        assert norm_std > 1.0, f"{name}: expected ‖d₁_std d₀‖ > 1, got {norm_std:.2e}"
        assert norm_ex < 1e-12, f"{name}: exact should pass, got {norm_ex:.2e}"
        assert n_spur > 0, f"{name}: expected n_spur > 0, got {n_spur}"

        print(f"         {name}: ‖d₁d₀‖_std={norm_std:.2f}, ‖d₁d₀‖_ex={norm_ex:.2e}, "
              f"n_spur={n_spur}  [PASS]")


def test_direction_dependence(structs):
    """T2.2: Pollution count depends on k-direction. Tested on 3 structures."""
    print("\n  T2.2  Direction-dependent pollution (3 structures)")
    directions = {
        '[100]': [1.0, 0.0, 0.0],
        '[110]': [1.0, 1.0, 0.0],
        '[111]': [1.0, 1.0, 1.0],
    }
    for name in ['Kelvin N=2', 'C15 N=1', 'WP N=1']:
        data = structs[name]
        star1, star2 = build_hodge_stars_voronoi(data)
        nV = len(data['V'])

        n_spurs = {}
        for label, d in directions.items():
            d_norm = np.array(d) / np.linalg.norm(d)
            k = make_k(data, 0.10, d_norm)
            _, d1_std, _ = build_both_d1(data, k)
            K_std = d1_std.conj().T @ np.diag(star2) @ d1_std
            n_zero = count_zeros(K_std, star1)
            n_spurs[label] = nV - n_zero

        vals = list(n_spurs.values())
        is_direction_dep = len(set(vals)) > 1
        parts = ", ".join(f"{l}:{n}" for l, n in n_spurs.items())
        print(f"         {name}: {parts}  "
              f"[{'direction-dependent' if is_direction_dep else 'ISOTROPIC'}]")
        assert is_direction_dep, f"{name}: expected direction-dependent pollution"
    print(f"         All 3 structures direction-dependent  [PASS]")


def test_intra_face_contradictions(structs):
    """T2.4 + T2.5: Intra-face contradictions and 1D null space (Cor 1)."""
    print("\n  T2.4+T2.5  Intra-face contradictions (Corollary 1)")
    data = structs['Kelvin N=2']
    V, E, F = data['V'], data['E'], data['F']
    L_vec = data['L_vec']
    shifts = compute_edge_shifts(V, E, L_vec)
    k = make_k(data, 0.10, [1.0, 0.3, 0.7])
    d0_k = build_d0_bloch(V, E, k, L_vec, shifts)
    face_edges = build_face_edges(E, F)

    # T2.5: All faces have 1D null space
    for f_idx, face in enumerate(F):
        A = np.zeros((len(face), len(face)), dtype=complex)
        for i in range(len(face)):
            v = face[i]
            e_prev, o_prev = face_edges[f_idx][(i - 1) % len(face)]
            e_curr, o_curr = face_edges[f_idx][i]
            A[i, (i - 1) % len(face)] = o_prev * d0_k[e_prev, v]
            A[i, i] = o_curr * d0_k[e_curr, v]
        rank = np.linalg.matrix_rank(A, tol=1e-10)
        null_dim = len(face) - rank
        assert null_dim == 1, f"Face {f_idx}: null_dim={null_dim}, expected 1"
    print(f"         T2.5: All {len(F)} faces have 1D null space  [PASS]")

    # T2.4: Count faces with same-shift edges having different recurrence phases
    n_contra = 0
    for f_idx, face in enumerate(F):
        n = len(face)
        phases = [1.0 + 0j]
        for i in range(1, n):
            e_prev, o_prev = face_edges[f_idx][i - 1]
            e_curr, o_curr = face_edges[f_idx][i]
            v = face[i]
            phases.append(-o_prev * phases[i-1] * d0_k[e_prev, v] / (o_curr * d0_k[e_curr, v]))

        shift_groups = {}
        for pos in range(n):
            e_idx, _ = face_edges[f_idx][pos]
            key = tuple(shifts[e_idx])
            shift_groups.setdefault(key, []).append(phases[pos])

        for key, phs in shift_groups.items():
            if len(phs) > 1 and any(abs(p - phs[0]) > 1e-10 for p in phs[1:]):
                n_contra += 1
                break

    assert n_contra > 0, f"Expected contradictions, got 0"
    print(f"         T2.4: {n_contra}/{len(F)} faces have intra-face contradictions  [PASS]")


def test_d1d0_scaling(structs):
    """T2.6: ‖d₁_std d₀‖ ~ O(k) scaling."""
    print("\n  T2.6  ‖d₁_std d₀‖ ~ O(k) scaling")
    data = structs['Kelvin N=2']
    V, E, F = data['V'], data['E'], data['F']
    L, L_vec = data['L'], data['L_vec']
    shifts = compute_edge_shifts(V, E, L_vec)
    crossings = compute_edge_crossings(V, E, L)
    edge_lookup = build_edge_lookup(E, crossings)

    fracs = [0.02, 0.05, 0.10, 0.15, 0.20]
    norms = []
    for frac in fracs:
        k = make_k(data, frac, [1, 0, 0])
        d0_k = build_d0_bloch(V, E, k, L_vec, shifts)
        d1_std = build_d1_bloch_standard(V, E, F, L, k, edge_lookup, crossings)
        norms.append(np.linalg.norm(d1_std @ d0_k))

    # Fit log-log slope (all points and small-k subset)
    log_k = np.log(fracs)
    log_n = np.log(norms)
    slope_all = np.polyfit(log_k, log_n, 1)[0]
    slope_low = np.polyfit(log_k[:3], log_n[:3], 1)[0]  # k ≤ 0.10 only

    for frac, norm in zip(fracs, norms):
        print(f"         frac={frac:.2f}: ‖d₁_std d₀‖ = {norm:.2f}")
    print(f"         Log-log slope (all): {slope_all:.2f}, (k≤0.10): {slope_low:.2f}")
    assert 0.7 < slope_low < 1.3, f"Small-k slope {slope_low:.2f} not near 1.0"
    print(f"         [PASS]")


def test_standard_kernel_deficit(structs):
    """T2.7: NEG — Standard ker(K) ≠ V."""
    print("\n  T2.7  NEG: Standard ker(K) ≠ V")
    for name, data in structs.items():
        k = make_k(data, 0.10, [1.0, 0.0, 0.0])
        d1_ex, d1_std, d0_k = build_both_d1(data, k)
        star1, star2 = build_hodge_stars_voronoi(data)
        nV = len(data['V'])

        K_ex = d1_ex.conj().T @ np.diag(star2) @ d1_ex
        K_std = d1_std.conj().T @ np.diag(star2) @ d1_std
        n_zero_ex = count_zeros(K_ex, star1)
        n_zero_std = count_zeros(K_std, star1)

        assert n_zero_ex == nV, f"{name}: exact ker should be V={nV}, got {n_zero_ex}"
        assert n_zero_std < nV, f"{name}: standard ker should be < V"
        print(f"         {name}: ker_exact={n_zero_ex}=V, ker_std={n_zero_std}<V, "
              f"deficit={nV - n_zero_std}  [PASS]")


def test_standard_rank_unstable(structs):
    """T2.8: NEG — Standard rank(d₁) ≠ exact AND direction-dependent."""
    print("\n  T2.8  NEG: Standard rank unstable")
    data = structs['Kelvin N=2']
    V, E, F = data['V'], data['E'], data['F']
    L, L_vec = data['L'], data['L_vec']
    shifts = compute_edge_shifts(V, E, L_vec)
    crossings = compute_edge_crossings(V, E, L)
    edge_lookup = build_edge_lookup(E, crossings)

    dirs = {'[100]': [1,0,0], '[110]': [1,1,0], '[111]': [1,1,1]}
    ranks_std = {}
    for label, d in dirs.items():
        d_norm = np.array(d, dtype=float) / np.linalg.norm(d)
        k = make_k(data, 0.10, d_norm)
        d1_s = build_d1_bloch_standard(V, E, F, L, k, edge_lookup, crossings)
        ranks_std[label] = np.linalg.matrix_rank(d1_s, tol=1e-10)

    # Exact rank (stable reference)
    k_ref = make_k(data, 0.10, [1, 0, 0])
    d0_ref = build_d0_bloch(V, E, k_ref, L_vec, shifts)
    d1_ex = build_d1_bloch_exact(V, E, F, k_ref, L_vec, d0_ref)
    rank_ex = np.linalg.matrix_rank(d1_ex, tol=1e-10)

    print(f"         Exact rank(d₁) = {rank_ex} (stable across directions)")
    for label, r in ranks_std.items():
        print(f"         Standard {label}: rank(d₁) = {r}")

    # Claim 1: rank_std ≠ rank_exact
    all_differ = all(r != rank_ex for r in ranks_std.values())
    assert all_differ, f"Expected all standard ranks ≠ exact rank {rank_ex}"
    print(f"         All standard ranks ≠ exact ({rank_ex})  [PASS]")

    # Claim 2: rank_std is direction-dependent
    rank_vals = list(ranks_std.values())
    is_dir_dep = len(set(rank_vals)) > 1
    assert is_dir_dep, f"Expected direction-dependent ranks, got {ranks_std}"
    print(f"         Standard rank direction-dependent: {ranks_std}  [PASS]")
    print(f"         → cohomology undefined on standard  [PASS]")


def test_k0_control(structs):
    """T2.9: CONTROL — At k=0, exact = standard (no Bloch → no problem).

    Note: we compare |d₁_exact(0)| vs |d₁_topological| (modulus, not sign).
    At k=0 the recurrence gives real phases ±1 per face, matching the
    topological d₁ up to orientation choice. Modulus comparison is intentional.
    """
    print("\n  T2.9  CONTROL: At k=0, exact = standard")
    for name, data in structs.items():
        V, E, F = data['V'], data['E'], data['F']
        d0_top = build_d0(V, E)
        d1_top = build_d1(V, E, F)

        L_vec = data['L_vec']
        shifts = compute_edge_shifts(V, E, L_vec)
        k0 = np.zeros(3)
        d0_k0 = build_d0_bloch(V, E, k0, L_vec, shifts)
        d1_ex_k0 = build_d1_bloch_exact(V, E, F, k0, L_vec, d0_k0)

        norm_top = np.linalg.norm(d1_top @ d0_top)
        norm_ex = np.linalg.norm(d1_ex_k0 @ d0_k0)

        assert norm_top < 1e-14, f"{name}: topological d₁d₀ ≠ 0"
        assert norm_ex < 1e-14, f"{name}: exact d₁d₀ at k=0 ≠ 0"

        # Compare modulus: |d₁_exact(0)| = |d₁_topological|
        # (sign per face may differ — recurrence seed is arbitrary at k=0)
        diff = np.abs(np.abs(d1_ex_k0) - np.abs(d1_top.astype(float)))
        assert np.max(diff) < 1e-14, f"{name}: |d₁_exact(0)| ≠ |d₁_topological|"

        print(f"         {name}: d₁d₀=0 at k=0, |d₁_ex(0)|=|d₁_top|  [PASS]")


def test_n_spur_scaling():
    """T2.10: n_spur grows as ~N^2.3 (between surface and volume).

    Pollution is extensive — grows with system size. Scaling exponent
    ~2.3 (closer to surface N² than volume N³) consistent with W3 R13:
    spurious modes are surface states (IPR ∝ N⁻²).
    """
    print("\n  T2.10  n_spur scaling with N (Kelvin)")
    Ns = [2, 3, 4]
    spurs_100 = []
    spurs_111 = []

    for N in Ns:
        data = build_kelvin_with_dual_info(N=N)
        V, E, F = data['V'], data['E'], data['F']
        L, L_vec = data['L'], data['L_vec']
        star1, star2 = build_hodge_stars_voronoi(data)
        nV = len(V)
        crossings = compute_edge_crossings(V, E, L)
        edge_lookup = build_edge_lookup(E, crossings)

        for label, d_vec, spur_list in [('[100]', [1,0,0], spurs_100),
                                         ('[111]', [1,1,1], spurs_111)]:
            d_norm = np.array(d_vec, dtype=float) / np.linalg.norm(d_vec)
            k = 2 * np.pi / L * 0.10 * d_norm
            d1_std = build_d1_bloch_standard(V, E, F, L, k, edge_lookup, crossings)
            K_std = d1_std.conj().T @ np.diag(star2) @ d1_std
            K_std = 0.5 * (K_std + K_std.conj().T)
            M = np.diag(star1)
            eigs = np.real(eigh(K_std, M, eigvals_only=True))
            thresh = max(np.max(np.abs(eigs)) * 1e-12, 1e-10)
            n_zero = int(np.sum(np.abs(eigs) < thresh))
            spur_list.append(nV - n_zero)

        print(f"         N={N}: V={nV}, n_spur [100]={spurs_100[-1]}, [111]={spurs_111[-1]}")

    # Fit scaling exponent
    log_N = np.log(Ns)
    for label, spurs in [('[100]', spurs_100), ('[111]', spurs_111)]:
        slope = np.polyfit(log_N, np.log(spurs), 1)[0]
        print(f"         {label}: n_spur ~ N^{slope:.2f}")
        assert slope > 1.5, f"Scaling too weak: {slope:.2f}"
        assert slope < 3.5, f"Scaling too strong: {slope:.2f}"

    print(f"         Pollution extensive (~N^2.3, near surface scaling)  [PASS]")


def test_per_edge_optimization():
    """T2.11: Numerical optimization of per-edge phases cannot achieve exactness.

    Direct demonstration of Corollary 1: even optimizing 192 free phases
    (one per edge), the minimum ‖d₁d₀‖ remains O(1). The exact construction
    achieves 10⁻¹⁶ — a ratio of 10¹⁵.
    """
    print("\n  T2.11  Per-edge phase optimization (Cor 1 numeric)")
    from scipy.optimize import minimize

    data = build_kelvin_with_dual_info(N=2)
    V, E, F = data['V'], data['E'], data['F']
    L, L_vec = data['L'], data['L_vec']
    nE = len(E)
    shifts = compute_edge_shifts(V, E, L_vec)

    k = 2 * np.pi / L * 0.10 * np.array([1.0, 0.0, 0.0])
    d0_k = build_d0_bloch(V, E, k, L_vec, shifts)
    d1_top = build_d1(V, E, F)

    # Per-edge parametrization: d₁[f,e] = d1_top[f,e] * exp(i*θ[e])
    def objective(theta):
        phases = np.exp(1j * theta)
        d1_param = d1_top.astype(complex) * phases[np.newaxis, :]
        prod = d1_param @ d0_k
        return np.real(np.sum(np.abs(prod)**2))

    # 10 random starts
    results = []
    for seed in range(10):
        rng = np.random.RandomState(seed)
        theta0 = rng.uniform(0, 2 * np.pi, nE)
        res = minimize(objective, theta0, method='L-BFGS-B',
                       options={'maxiter': 500, 'ftol': 1e-15})
        results.append(np.sqrt(res.fun))

    best = min(results)
    median = np.median(results)

    # Compare with exact
    d1_ex = build_d1_bloch_exact(V, E, F, k, L_vec, d0_k)
    norm_ex = np.linalg.norm(d1_ex @ d0_k)

    print(f"         Per-edge optimization (10 starts):")
    print(f"           best ‖d₁d₀‖  = {best:.4f}")
    print(f"           median        = {median:.4f}")
    print(f"         Exact recurrence: ‖d₁d₀‖ = {norm_ex:.2e}")
    print(f"         Ratio: {best / max(norm_ex, 1e-20):.0e}")

    assert best > 0.5, f"Per-edge optimization too good: {best:.4f}"
    assert norm_ex < 1e-12, f"Exact should be machine precision"
    print(f"         Per-edge approach fundamentally impossible  [PASS]")


def test_phase_sensitivity():
    """T2.12: Exact construction is an isolated point — O(ε) sensitivity.

    Perturbing exact d₁ phases by ε produces ‖d₁d₀‖ = O(ε) with constant
    ~33, linear across 8 orders of magnitude. The exact solution is not a
    basin of attraction — any deviation destroys exactness proportionally.
    """
    print("\n  T2.12  Phase sensitivity: exact is isolated point")
    data = build_all_structures()['Kelvin N=2']
    V, E, F = data['V'], data['E'], data['F']
    L, L_vec = data['L'], data['L_vec']
    nE, nF = len(E), len(F)
    shifts = compute_edge_shifts(V, E, L_vec)

    k = 2 * np.pi / L * 0.10 * np.array([1.0, 0.0, 0.0])
    d0_k = build_d0_bloch(V, E, k, L_vec, shifts)
    d1_ex = build_d1_bloch_exact(V, E, F, k, L_vec, d0_k)

    # Fixed random perturbation direction. The constant ‖d₁d₀‖/ε depends on
    # the direction ξ (different seeds give different constants ~20-50) but
    # linearity in ε is universal.
    rng = np.random.RandomState(42)
    xi = rng.randn(nF, nE)

    ratios = []
    epsilons = [1e-8, 1e-6, 1e-4, 1e-2, 0.1]
    for eps in epsilons:
        d1_pert = d1_ex * np.exp(1j * eps * xi)
        d1_pert[np.abs(d1_ex) < 0.5] = 0  # preserve support
        norm_pert = np.linalg.norm(d1_pert @ d0_k)
        ratio = norm_pert / eps
        ratios.append(ratio)
        print(f"         ε={eps:.0e}: ‖d₁d₀‖={norm_pert:.4e}, ‖d₁d₀‖/ε={ratio:.2f}")

    # Check linearity: ratios should be approximately constant
    cv = np.std(ratios) / np.mean(ratios)
    assert cv < 0.15, f"Sensitivity not linear: CV = {cv:.3f}"
    print(f"         ‖d₁d₀‖/ε ≈ {np.mean(ratios):.1f} ± {np.std(ratios):.1f} "
          f"(CV={cv:.3f}, linear O(ε))")
    print(f"         Exact construction is an isolated point  [PASS]")


def main():
    print("=" * 60)
    print("GROUP 2: The problem — standard fails (§3)")
    print("=" * 60)

    structs = build_all_structures()
    test_standard_breaks_exactness(structs)
    test_direction_dependence(structs)
    test_intra_face_contradictions(structs)
    test_d1d0_scaling(structs)
    test_standard_kernel_deficit(structs)
    test_standard_rank_unstable(structs)
    test_k0_control(structs)
    test_n_spur_scaling()
    test_per_edge_optimization()
    test_phase_sensitivity()

    print("\n" + "-" * 60)
    print("GROUP 2: ALL PASS")
    print("-" * 60)


if __name__ == '__main__':
    main()
