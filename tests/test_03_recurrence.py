"""
GROUP 3 — The fix: recurrence construction (§4).

Verify the construction and its algebraic properties.

Tests:
  T3.1   ‖d₁_exact(k) d₀(k)‖ < 10⁻¹⁴ on all structures, multiple k     [Thm 1]
  T3.2   Holonomy H_f = 1 on all faces, all structures                    [Lem 1]
  T3.3   K canonical under face permutations (6 types)                    [Prop 2]
  T3.4   K invariant under gauge transform e^{iθ(v)}                     [Prop 2]
  T3.5   n_zero = V at all k, all structures                             [Prop 3]
  T3.6   Rank chain: rank(d₀) = V, ker(d₁) = im(d₀), max grad% < 10⁻⁸  [Prop 3]
  T3.7   Exactness preserved under mesh perturbation ε = 0..20%          [Thm 1]
  T3.8   d₂d₁ = 0 (PARTIAL: d₂ Bloch not implemented, topological only)  [§4.6]
  T3.9   β(Γ) = (1,3,3,1), β₀=β₁=0 at k≠0 (PARTIAL: β₂,β₃ need d₂(k)) [Remark]
  T3.10  BZ boundary: exactness maintained, rank drop expected            [Nice]

Note: T3.3 reverse-invariance is trivial (K=d₁†⋆₂d₁ absorbs sign).
      T3.7 tests algebraic exactness only, not metric consistency.
      T3.8/T3.9 are PARTIAL — d₂ Bloch recurrence needed for completion.

Structures: Kelvin N=2, C15 N=1, WP N=1, SC N=3.

Expected output (verified 17 Mar 2026):

  ============================================================
  GROUP 3: The fix — recurrence construction (§4)
  ============================================================

    T3.1  Exactness: ‖d₁_exact d₀‖ < 10⁻¹⁴
           Kelvin N=2: max ‖d₁d₀‖ = 8.29e-16  (8 k-points)  [PASS]
           C15 N=1: max ‖d₁d₀‖ = 9.53e-16  (8 k-points)  [PASS]
           WP N=1: max ‖d₁d₀‖ = 6.93e-16  (8 k-points)  [PASS]
           SC N=3: max ‖d₁d₀‖ = 8.45e-16  (8 k-points)  [PASS]

    T3.2  Holonomy H_f = 1 (Lemma 1)
           Kelvin N=2: max |H_f - 1| = 2.22e-16  (112 faces)  [PASS]
           C15 N=1: max |H_f - 1| = 2.30e-16  (160 faces)  [PASS]
           WP N=1: max |H_f - 1| = 2.22e-16  (54 faces)  [PASS]
           SC N=3: max |H_f - 1| = 2.22e-16  (81 faces)  [PASS]

    T3.3  Canonical K under permutations (Prop 2)
           cyclic_1: ‖K_perm − K_ref‖ = 2.80e-15  [PASS]
           cyclic_2: ‖K_perm − K_ref‖ = 5.17e-15  [PASS]
           reverse: ‖K_perm − K_ref‖ = 4.74e-15  [PASS]

    T3.4  K invariant under gauge transform (Prop 2)
           5 random gauges: max |Δeig| < 10⁻¹²  [PASS]

    T3.5  ker(K) = V at all k (Prop 3)
           Kelvin N=2: n_zero = 96 = V at all 8 k-points  [PASS]
           C15 N=1: n_zero = 136 = V at all 8 k-points  [PASS]
           WP N=1: n_zero = 46 = V at all 8 k-points  [PASS]
           SC N=3: n_zero = 27 = V at all 8 k-points  [PASS]

    T3.6  Rank chain: im(d₀) = ker(d₁) = ker(K) (Prop 3)
           Kelvin N=2: rank(d₀)=96, ker(d₁)=96, ker(K)=96, max_grad=6.41e-14  [PASS]
           C15 N=1: rank(d₀)=136, ker(d₁)=136, ker(K)=136, max_grad=1.34e-13  [PASS]
           WP N=1: rank(d₀)=46, ker(d₁)=46, ker(K)=46, max_grad=5.63e-14  [PASS]
           SC N=3: rank(d₀)=27, ker(d₁)=27, ker(K)=27, max_grad=4.56e-14  [PASS]

    T3.7  Exactness under perturbation ε = 0..20% (algebraic, not metric)
           ε=0%: ‖d₁d₀‖ = 7.86e-16  [PASS]
           ε=5%: ‖d₁d₀‖ = 7.86e-16  [PASS]
           ε=10%: ‖d₁d₀‖ = 7.86e-16  [PASS]
           ε=15%: ‖d₁d₀‖ = 7.86e-16  [PASS]
           ε=20%: ‖d₁d₀‖ = 7.86e-16  [PASS]

    T3.8  d₂d₁ = 0 (PARTIAL — d₂ Bloch recurrence not implemented)
           SC N=3: ‖d₂_top d₁_top‖ = 0.00e+00, ‖d₂_top d₁_exact‖ = 3.71e+00
           d₂ also needs Bloch recurrence (same pattern as d₁)  [PARTIAL]

    T3.9  Bloch cohomology
           SC N=3: β(Γ) = (1, 3, 3, 1)  [PASS]
           SC N=3: β₀(k≠0)=0, β₁(k≠0)=0  (β₂,β₃ need d₂(k))  [PARTIAL]

    T3.10  BZ boundary behavior
           frac=0.500: ‖d₁d₀‖=6.36e-32, rank(d₀)=96, n_zero=96
           frac=0.900: ‖d₁d₀‖=4.82e-16, rank(d₀)=96, n_zero=96
           frac=0.999: ‖d₁d₀‖=2.56e-16, rank(d₀)=96, n_zero=96
           frac=1.000: ‖d₁d₀‖=2.55e-31, rank(d₀)=95, n_zero=98
           Exactness maintained at BZ boundary  [PASS]

  ------------------------------------------------------------
  GROUP 3: 8 PASS, 2 PARTIAL (0.8s)
  ------------------------------------------------------------
"""

import numpy as np
from scipy.linalg import eigh

from physics.hodge import (
    build_kelvin_with_dual_info,
    build_c15_with_dual_info,
    build_wp_with_dual_info,
    build_hodge_stars_voronoi,
)
from physics.gauge_bloch import compute_edge_shifts, build_d0_bloch, build_d1_bloch_exact
from core_math.builders.solids_periodic import build_sc_supercell_periodic
from core_math.operators.incidence import build_d0, build_d1
from physics.bloch_complex import (
    build_face_edge_map, build_d2_top_from_foam, build_d2_bloch_exact,
)


# ── Helpers ──────────────────────────────────────────────────────────────

def wrap_sc(N):
    a = 2.0
    L = a * N
    verts, edges, faces, cells = build_sc_supercell_periodic(N=N)
    return {'V': verts, 'E': edges, 'F': faces, 'C': cells,
            'L': L, 'L_vec': np.array([L, L, L])}


def build_structures():
    return {
        'Kelvin N=2': build_kelvin_with_dual_info(N=2),
        'C15 N=1':    build_c15_with_dual_info(N=1),
        'WP N=1':     build_wp_with_dual_info(N=1),
        'SC N=3':     wrap_sc(N=3),
    }


def make_k(data, frac, direction):
    d = np.array(direction, dtype=float)
    d /= np.linalg.norm(d)
    return 2 * np.pi / data['L'] * frac * d


def build_K_M(d1_k, star1, star2):
    K = d1_k.conj().T @ np.diag(star2) @ d1_k
    K = 0.5 * (K + K.conj().T)
    M = np.diag(star1)
    return K, M


def count_zeros(eigs):
    thresh = max(np.max(np.abs(eigs)) * 1e-12, 1e-10)
    return int(np.sum(np.abs(eigs) < thresh))


def get_star(data):
    """Get Hodge stars. SC uses uniform stars."""
    if 'cell_centers' in data:
        return build_hodge_stars_voronoi(data)
    # SC: uniform edge length a=2, face area a², volume a³
    nE, nF = len(data['E']), len(data['F'])
    a = 2.0
    return np.full(nE, a), np.full(nF, 1.0 / a)


# ── Tests ────────────────────────────────────────────────────────────────

def test_exactness_all_structures(structs):
    """T3.1: ‖d₁_exact d₀‖ < 10⁻¹⁴ on all structures, multiple k-points."""
    print("\n  T3.1  Exactness: ‖d₁_exact d₀‖ < 10⁻¹⁴")
    for name, data in structs.items():
        V, E, F = data['V'], data['E'], data['F']
        L_vec = data['L_vec']
        shifts = compute_edge_shifts(V, E, L_vec)

        max_norm = 0
        for frac in [0.02, 0.05, 0.10, 0.20]:
            for d in [[1, 0, 0], [1, 1, 1]]:
                k = make_k(data, frac, d)
                d0_k = build_d0_bloch(V, E, k, L_vec, shifts)
                d1_k = build_d1_bloch_exact(V, E, F, k, L_vec, d0_k)
                norm = np.linalg.norm(d1_k @ d0_k)
                max_norm = max(max_norm, norm)

        assert max_norm < 1e-14, f"{name}: max ‖d₁d₀‖ = {max_norm:.2e}"
        print(f"         {name}: max ‖d₁d₀‖ = {max_norm:.2e}  (8 k-points)  [PASS]")


def test_holonomy(structs):
    """T3.2: Holonomy H_f = 1 on all faces, all structures.

    H_f = ∏ᵢ (-σᵢ₋₁ d₀[eᵢ₋₁,vᵢ]) / (σᵢ d₀[eᵢ,vᵢ])
    If H_f = 1, the recurrence closes (last vertex equation automatic).
    """
    print("\n  T3.2  Holonomy H_f = 1 (Lemma 1)")
    for name, data in structs.items():
        V, E, F = data['V'], data['E'], data['F']
        L_vec = data['L_vec']
        shifts = compute_edge_shifts(V, E, L_vec)

        edge_map = {}
        for idx, (i, j) in enumerate(E):
            edge_map[(i, j)] = (idx, +1)
            edge_map[(j, i)] = (idx, -1)

        k = make_k(data, 0.10, [1.0, 0.3, 0.7])
        d0_k = build_d0_bloch(V, E, k, L_vec, shifts)

        max_err = 0
        for face in F:
            n = len(face)
            # Build edge info for this face
            edges_info = []
            for v_pos in range(n):
                i, j = face[v_pos], face[(v_pos + 1) % n]
                e_idx, orient = edge_map[(i, j)]
                edges_info.append((e_idx, orient))

            # Compute holonomy: product of ratios around face
            H = 1.0 + 0j
            for i in range(n):
                v = face[(i + 1) % n]
                e_curr, o_curr = edges_info[i]
                e_next, o_next = edges_info[(i + 1) % n]
                # Ratio at vertex v: -σ_curr * d0[e_curr, v] / (σ_next * d0[e_next, v])
                num = -o_curr * d0_k[e_curr, v]
                den = o_next * d0_k[e_next, v]
                H *= num / den

            max_err = max(max_err, abs(H - 1.0))

        assert max_err < 1e-13, f"{name}: max |H_f - 1| = {max_err:.2e}"
        print(f"         {name}: max |H_f - 1| = {max_err:.2e}  ({len(F)} faces)  [PASS]")


def test_canonical_K_permutations(structs):
    """T3.3: K is canonical under face vertex permutations."""
    print("\n  T3.3  Canonical K under permutations (Prop 2)")
    data = structs['Kelvin N=2']
    V, E, F = data['V'], data['E'], data['F']
    L_vec = data['L_vec']
    star1, star2 = get_star(data)
    shifts = compute_edge_shifts(V, E, L_vec)
    k = make_k(data, 0.10, [1.0, 0.0, 0.0])
    d0_k = build_d0_bloch(V, E, k, L_vec, shifts)

    # Reference K
    d1_ref = build_d1_bloch_exact(V, E, F, k, L_vec, d0_k)
    K_ref = d1_ref.conj().T @ np.diag(star2) @ d1_ref

    # Permute faces: cyclic shift and reverse (both valid cyclic orderings).
    # Random permutation is INVALID — breaks consecutive vertex adjacency.
    # Note: reverse changes face orientation (sign of d₁ row), but K = d₁†⋆₂d₁
    # absorbs the sign (|λ_f|²=1). So reverse-invariance is a consequence of
    # K's quadratic structure, not of the recurrence itself.
    permutations = {
        'cyclic_1': lambda f: f[1:] + f[:1],
        'cyclic_2': lambda f: f[2:] + f[:2],
        'reverse':  lambda f: f[::-1],
    }

    for perm_name, perm_fn in permutations.items():
        F_perm = [perm_fn(list(face)) for face in F]
        d1_perm = build_d1_bloch_exact(V, E, F_perm, k, L_vec, d0_k)
        K_perm = d1_perm.conj().T @ np.diag(star2) @ d1_perm
        diff = np.linalg.norm(K_perm - K_ref)
        assert diff < 1e-12, f"{perm_name}: ‖ΔK‖ = {diff:.2e}"
        print(f"         {perm_name}: ‖K_perm − K_ref‖ = {diff:.2e}  [PASS]")


def test_gauge_invariance(structs):
    """T3.4: K invariant under vertex gauge transform e^{iθ(v)}."""
    print("\n  T3.4  K invariant under gauge transform (Prop 2)")
    data = structs['Kelvin N=2']
    V, E, F = data['V'], data['E'], data['F']
    L_vec = data['L_vec']
    star1, star2 = get_star(data)
    shifts = compute_edge_shifts(V, E, L_vec)
    k = make_k(data, 0.10, [1.0, 0.0, 0.0])
    d0_k = build_d0_bloch(V, E, k, L_vec, shifts)

    d1_ref = build_d1_bloch_exact(V, E, F, k, L_vec, d0_k)
    K_ref, M = build_K_M(d1_ref, star1, star2)
    eigs_ref = np.sort(np.real(eigh(K_ref, M, eigvals_only=True)))

    rng = np.random.RandomState(123)
    for seed in range(5):
        theta = rng.uniform(0, 2 * np.pi, len(V))
        # Gauge transform: d0 → D†d0D where D = diag(e^{iθ})
        D = np.diag(np.exp(1j * theta))
        D_E = np.diag(np.exp(1j * np.array([
            theta[E[e][1]] - theta[E[e][0]] for e in range(len(E))
        ])))
        # d0_k transforms, d1 adapts, K should be same
        d0_gauged = d0_k @ D
        d1_gauged = build_d1_bloch_exact(V, E, F, k, L_vec, d0_gauged)
        K_gauged, _ = build_K_M(d1_gauged, star1, star2)
        eigs_gauged = np.sort(np.real(eigh(K_gauged, M, eigvals_only=True)))

        max_diff = np.max(np.abs(eigs_gauged - eigs_ref))
        assert max_diff < 1e-12, f"seed {seed}: max |Δeig| = {max_diff:.2e}"

    print(f"         5 random gauges: max |Δeig| < 10⁻¹²  [PASS]")


def test_kernel_equals_V(structs):
    """T3.5: n_zero = V at all k, all structures."""
    print("\n  T3.5  ker(K) = V at all k (Prop 3)")
    for name, data in structs.items():
        V, E, F = data['V'], data['E'], data['F']
        L_vec = data['L_vec']
        star1, star2 = get_star(data)
        shifts = compute_edge_shifts(V, E, L_vec)
        nV = len(V)

        all_ok = True
        tested = 0
        for frac in [0.02, 0.05, 0.10, 0.20]:
            for d in [[1, 0, 0], [1, 1, 1]]:
                k = make_k(data, frac, d)
                d0_k = build_d0_bloch(V, E, k, L_vec, shifts)
                d1_k = build_d1_bloch_exact(V, E, F, k, L_vec, d0_k)
                K, M = build_K_M(d1_k, star1, star2)
                eigs = np.real(eigh(K, M, eigvals_only=True))
                n_zero = count_zeros(eigs)
                if n_zero != nV:
                    all_ok = False
                tested += 1

        assert all_ok, f"{name}: ker(K) ≠ V at some k"
        print(f"         {name}: n_zero = {nV} = V at all {tested} k-points  [PASS]")


def test_rank_chain(structs):
    """T3.6: rank(d₀) = V, ker(d₁) = im(d₀), physical modes gradient-free."""
    print("\n  T3.6  Rank chain: im(d₀) = ker(d₁) = ker(K) (Prop 3)")
    for name, data in structs.items():
        V, E, F = data['V'], data['E'], data['F']
        L_vec = data['L_vec']
        star1, star2 = get_star(data)
        shifts = compute_edge_shifts(V, E, L_vec)
        nV = len(V)
        M = np.diag(star1)

        k = make_k(data, 0.10, [1, 0, 0])
        d0_k = build_d0_bloch(V, E, k, L_vec, shifts)
        d1_k = build_d1_bloch_exact(V, E, F, k, L_vec, d0_k)

        rank_d0 = np.linalg.matrix_rank(d0_k)
        dim_ker_d1 = len(E) - np.linalg.matrix_rank(d1_k)

        K, _ = build_K_M(d1_k, star1, star2)
        eigs, vecs = eigh(K, M)
        eigs = np.real(eigs)
        idx = np.argsort(eigs)
        eigs, vecs = eigs[idx], vecs[:, idx]
        thresh = max(np.max(np.abs(eigs)) * 1e-12, 1e-10)
        dim_ker_K = int(np.sum(np.abs(eigs) < thresh))

        assert rank_d0 == nV, f"{name}: rank(d₀) = {rank_d0} ≠ V={nV}"
        assert dim_ker_d1 == nV, f"{name}: dim ker(d₁) = {dim_ker_d1} ≠ V={nV}"
        assert dim_ker_K == nV, f"{name}: dim ker(K) = {dim_ker_K} ≠ V={nV}"

        # Max gradient fraction of physical eigenvectors
        phys = vecs[:, np.abs(eigs) > thresh]
        if phys.shape[1] > 0:
            Md0 = M @ d0_k
            Pg = d0_k @ np.linalg.solve(d0_k.conj().T @ Md0, Md0.conj().T)
            max_g = 0
            for j in range(phys.shape[1]):
                v = phys[:, j]
                Pv = Pg @ v
                g = np.sqrt(abs(np.real(Pv.conj() @ M @ Pv) /
                                np.real(v.conj() @ M @ v)))
                max_g = max(max_g, g)
        else:
            max_g = 0

        assert max_g < 1e-8, f"{name}: max gradient fraction = {max_g:.2e}"
        print(f"         {name}: rank(d₀)={rank_d0}, ker(d₁)={dim_ker_d1}, "
              f"ker(K)={dim_ker_K}, max_grad={max_g:.2e}  [PASS]")


def test_perturbation_stability(structs):
    """T3.7: Exactness preserved under mesh perturbation ε = 0..50%.

    Note: this tests ALGEBRAIC exactness only (d₁d₀=0). Hodge stars are not
    rebuilt from perturbed geometry — the point is that exactness is a purely
    topological property independent of metric (Theorem 1 depends only on
    incidence structure and Bloch phases, not on edge lengths or face areas).
    """
    print("\n  T3.7  Exactness under perturbation ε = 0..50% (algebraic, not metric)")
    data = structs['Kelvin N=2']
    V, E, F = data['V'], data['E'], data['F']
    L_vec = data['L_vec']

    rng = np.random.RandomState(42)
    k = make_k(data, 0.10, [1, 0, 0])

    # Mean edge length for scaling
    edge_lengths = [np.linalg.norm(
        V[e[1]] - V[e[0]] - np.round((V[e[1]] - V[e[0]]) / L_vec) * L_vec
    ) for e in E]
    mean_len = np.mean(edge_lengths)

    for eps_frac in [0.0, 0.05, 0.10, 0.15, 0.20, 0.50]:
        eps = eps_frac * mean_len
        V_pert = V + rng.randn(*V.shape) * eps
        # Rebuild with perturbed vertices (same topology)
        data_pert = dict(data)
        data_pert['V'] = V_pert

        shifts_pert = compute_edge_shifts(V_pert, E, L_vec)
        d0_k = build_d0_bloch(V_pert, E, k, L_vec, shifts_pert)
        d1_k = build_d1_bloch_exact(V_pert, E, F, k, L_vec, d0_k)
        norm = np.linalg.norm(d1_k @ d0_k)
        assert norm < 1e-12, f"ε={eps_frac}: ‖d₁d₀‖ = {norm:.2e}"
        print(f"         ε={eps_frac:.0%}: ‖d₁d₀‖ = {norm:.2e}  [PASS]")


def test_d2_exactness(structs):
    """T3.8: d₂(k)d₁(k) = 0 — recurrence extends to level 2→3.

    Same recurrence pattern as d₁: derive face phases from d₁(k) entries
    via BFS on cell boundaries. Verified on SC (has cell data) and foams
    (build cell-face incidence from face_to_cells).
    """
    print("\n  T3.8  d₂(k)d₁(k) = 0 (full complex, level 2→3)")

    for name, data in structs.items():
        V, E, F = data['V'], data['E'], data['F']
        L_vec = data['L_vec']
        shifts = compute_edge_shifts(V, E, L_vec)
        nF = len(F)
        face_edges = build_face_edge_map(F, E)

        # Determine cell-face incidence
        if 'C' in data:
            # SC: cell data directly available
            C = data['C']
            nC = len(C)
            cfi = C  # already list of [(f_idx, orient), ...]
            # Build d₂_top for comparison
            d2_top = np.zeros((nC, nF), dtype=float)
            for c_idx, cell in enumerate(C):
                for f_idx, orient in cell:
                    d2_top[c_idx, f_idx] = orient
        elif 'face_to_cells' in data:
            # Foam: build from face_to_cells
            d2_top, cfi = build_d2_top_from_foam(data)
            nC = d2_top.shape[0]
        else:
            print(f"         {name}: no cell data, skip")
            continue

        # Verify d₂_top · d₁(k=0) = 0
        k0 = np.zeros(3)
        d0_k0 = build_d0_bloch(V, E, k0, L_vec, shifts)
        d1_k0 = build_d1_bloch_exact(V, E, F, k0, L_vec, d0_k0)
        norm_k0 = np.linalg.norm(d2_top @ d1_k0)
        assert norm_k0 < 1e-10, f"{name}: d₂_top·d₁(k=0) ≠ 0: {norm_k0:.2e}"

        # At k≠0: d₂_top · d₁_exact ≠ 0 (standard fails)
        k = make_k(data, 0.10, [1, 0, 0])
        d0_k = build_d0_bloch(V, E, k, L_vec, shifts)
        d1_k = build_d1_bloch_exact(V, E, F, k, L_vec, d0_k)
        norm_top_bloch = np.linalg.norm(d2_top.astype(complex) @ d1_k)

        # d₂ exact via recurrence
        d2_ex = build_d2_bloch_exact(cfi, face_edges, d1_k, nC, nF)
        norm_exact = np.linalg.norm(d2_ex @ d1_k)

        # Unimodular check
        unimod = np.allclose(np.abs(d2_ex[np.abs(d2_ex) > 0.5]), 1.0, atol=1e-10)

        print(f"         {name} (C={nC}): "
              f"‖d₂_exact·d₁‖={norm_exact:.2e}, "
              f"‖d₂_top·d₁‖={norm_top_bloch:.2e}, "
              f"unimod={unimod}")

        assert norm_exact < 1e-12, \
            f"{name}: d₂_exact·d₁_exact ≠ 0: {norm_exact:.2e}"
        assert norm_top_bloch > 0.5, \
            f"{name}: d₂_top should fail at k≠0: {norm_top_bloch:.2e}"
        assert unimod, f"{name}: d₂_exact entries should be unimodular"

    print("         Full complex: d₂(k)d₁(k) = 0 via recurrence  [PASS]")


def test_cohomology(structs):
    """T3.9: β(Γ) = (1,3,3,1), β(k≠0) = (0,0,0,0)."""
    print("\n  T3.9  Bloch cohomology")
    for name, data in structs.items():
        V, E, F = data['V'], data['E'], data['F']
        L_vec = data['L_vec']
        shifts = compute_edge_shifts(V, E, L_vec)
        nV, nE, nF = len(V), len(E), len(F)
        face_edges = build_face_edge_map(F, E)

        # Determine cell-face incidence
        if 'C' in data:
            C = data['C']
            nC = len(C)
            cfi = C
            d2_top = np.zeros((nC, nF), dtype=float)
            for c_idx, cell in enumerate(C):
                for f_idx, orient in cell:
                    d2_top[c_idx, f_idx] = orient
        elif 'face_to_cells' in data:
            d2_top, cfi = build_d2_top_from_foam(data)
            nC = d2_top.shape[0]
        else:
            print(f"         {name}: no cell data, skip")
            continue

        # ── β at Γ ──
        k0 = np.zeros(3)
        d0_k0 = build_d0_bloch(V, E, k0, L_vec, shifts)
        d1_k0 = build_d1_bloch_exact(V, E, F, k0, L_vec, d0_k0)

        r0 = np.linalg.matrix_rank(d0_k0, tol=1e-10)
        r1 = np.linalg.matrix_rank(d1_k0, tol=1e-10)
        r2 = np.linalg.matrix_rank(d2_top, tol=1e-10)

        b0 = nV - r0
        b1 = nE - r0 - r1
        b2 = nF - r1 - r2
        b3 = nC - r2

        betti_gamma = (int(b0), int(b1), int(b2), int(b3))
        assert betti_gamma == (1, 3, 3, 1), \
            f"{name}: β(Γ) = {betti_gamma}, expected (1,3,3,1)"
        print(f"         {name}: β(Γ) = {betti_gamma}  [PASS]")

        # ── β at generic k ≠ 0: full verification with d₂(k) ──
        k_gen = make_k(data, 0.10, [1, 0.3, 0.7])
        d0_kg = build_d0_bloch(V, E, k_gen, L_vec, shifts)
        d1_kg = build_d1_bloch_exact(V, E, F, k_gen, L_vec, d0_kg)
        d2_kg = build_d2_bloch_exact(cfi, face_edges, d1_kg, nC, nF)

        r0g = np.linalg.matrix_rank(d0_kg, tol=1e-10)
        r1g = np.linalg.matrix_rank(d1_kg, tol=1e-10)
        r2g = np.linalg.matrix_rank(d2_kg, tol=1e-10)

        b0g = int(nV - r0g)
        b1g = int(nE - r0g - r1g)
        b2g = int(nF - r1g - r2g)
        b3g = int(nC - r2g)

        betti_k = (b0g, b1g, b2g, b3g)
        assert betti_k == (0, 0, 0, 0), \
            f"{name}: β(k≠0) = {betti_k}, expected (0,0,0,0)"
        print(f"         {name}: β(k≠0) = {betti_k}  [PASS]")


def test_cohomology_trim(structs):
    """T3.12: Cohomology at TRIM points matches Künneth exactly.

    At TRIM points (X, M, R), some monodromies = -1. Since -1 ≠ 1,
    the corresponding S¹ factor is acyclic → β = (0,0,0,0).
    At BZ boundary k = (2π/L, 0, 0), all monodromies = 1 (equivalent
    to Γ) → β₁ = 3, n_zero > V.
    """
    print("\n  T3.12  Cohomology at TRIM points (Künneth verification)")
    data = structs['Kelvin N=2']
    V, E, F = data['V'], data['E'], data['F']
    L, L_vec = data['L'], data['L_vec']
    nV, nE = len(V), len(E)
    star1, star2 = get_star(data)
    shifts = compute_edge_shifts(V, E, L_vec)

    trim_points = {
        'X=(π/L,0,0)':     np.array([0.5, 0.0, 0.0]),
        'M=(π/L,π/L,0)':   np.array([0.5, 0.5, 0.0]),
        'R=(π/L,π/L,π/L)': np.array([0.5, 0.5, 0.5]),
        'BZ_x=(2π/L,0,0)': np.array([1.0, 0.0, 0.0]),
    }

    for label, frac_vec in trim_points.items():
        k = 2 * np.pi / L * frac_vec
        monodromies = [np.exp(1j * k[j] * L_vec[j]) for j in range(3)]
        m_trivial = sum(1 for lam in monodromies if abs(lam - 1) < 1e-10)

        d0_k = build_d0_bloch(V, E, k, L_vec, shifts)
        d1_k = build_d1_bloch_exact(V, E, F, k, L_vec, d0_k)

        rank_d0 = np.linalg.matrix_rank(d0_k, tol=1e-10)
        dim_ker_d1 = nE - np.linalg.matrix_rank(d1_k, tol=1e-10)
        beta1 = dim_ker_d1 - rank_d0

        K, M = build_K_M(d1_k, star1, star2)
        eigs = np.real(eigh(K, M, eigvals_only=True))
        n_zero = count_zeros(eigs)

        # Künneth on T³ = S¹×S¹×S¹: H^p(T³,L_k) = ⊗_j H^*(S¹_j, L_j).
        # If ANY factor has nontrivial monodromy (λ_j ≠ 1), that S¹ is
        # acyclic (H⁰=H¹=0), so the tensor product vanishes entirely.
        # β₁ = 3 only when ALL monodromies = 1 (equivalent to Γ).
        # Cross-validates T3.10: n_zero=98 at BZ_x is this Γ-equivalence.
        all_trivial = (m_trivial == 3)
        beta1_pred = 3 if all_trivial else 0

        ok = beta1 == beta1_pred
        print(f"         {label}: m={m_trivial}, β₁={beta1} (pred={beta1_pred}), "
              f"n_zero={n_zero}  [{'PASS' if ok else 'FAIL'}]")
        assert ok, f"{label}: β₁={beta1} ≠ predicted {beta1_pred}"

    print(f"         All TRIM points match Künneth exactly  [PASS]")


def test_iff_completeness(structs):
    """T3.11: Solution space of d₁d₀=0 is EXACTLY 1D per face (iff).

    The recurrence family is the ONLY set of local operators satisfying
    exactness. No other d₁ with support on ∂f exists outside this family.

    Method: build the full constraint matrix A (d₁d₀=0 linearized over
    all nonzero entries of d₁) and verify that nullity(A) = nF.
    Since each face contributes exactly 1 free parameter (the seed φ₀),
    this confirms the solution space is 1D per face — hence "iff."
    """
    print("\n  T3.11  iff: solution space = exactly 1D per face")
    from core_math.operators.incidence import build_d1 as build_d1_top

    for name, data in structs.items():
        V, E, F = data['V'], data['E'], data['F']
        L_vec = data['L_vec']
        nV, nE, nF = len(V), len(E), len(F)
        shifts = compute_edge_shifts(V, E, L_vec)
        d1_top = build_d1_top(V, E, F)

        k = make_k(data, 0.10, [1.0, 0.3, 0.7])
        d0_k = build_d0_bloch(V, E, k, L_vec, shifts)

        # Build variable map: one complex variable per nonzero entry of d₁
        var_map = {}
        idx = 0
        for f in range(nF):
            for e in range(nE):
                if abs(d1_top[f, e]) > 0:
                    var_map[(f, e)] = idx
                    idx += 1
        total_vars = idx

        # Build constraint matrix A: for each (f, v), equation
        # Σ_e d₁[f,e] · d₀(k)[e,v] = 0 with e ∈ ∂f, v ∈ e
        # Iterate only over edges in ∂f (not all nE edges)
        rows = []
        for f in range(nF):
            face = F[f]
            face_edges = [e for e in range(nE) if (f, e) in var_map]
            for v in face:
                row = np.zeros(total_vars, dtype=complex)
                for e in face_edges:
                    if abs(d0_k[e, v]) > 0:
                        row[var_map[(f, e)]] = d0_k[e, v]
                rows.append(row)

        A = np.array(rows)
        rank_A = np.linalg.matrix_rank(A, tol=1e-10)
        nullity = total_vars - rank_A

        assert nullity == nF, \
            f"{name}: nullity = {nullity} ≠ nF = {nF}"
        print(f"         {name}: vars={total_vars}, rank={rank_A}, "
              f"nullity={nullity} = nF={nF}  [PASS]")

    print(f"         Recurrence family is the ONLY solution (iff)  [PASS]")


def test_bz_boundary(structs):
    """T3.10: BZ boundary — exactness maintained, rank drop expected."""
    print("\n  T3.10  BZ boundary behavior")
    data = structs['Kelvin N=2']
    V, E, F = data['V'], data['E'], data['F']
    L_vec = data['L_vec']
    star1, star2 = get_star(data)
    shifts = compute_edge_shifts(V, E, L_vec)
    nV = len(V)

    for frac in [0.50, 0.90, 0.999, 1.0]:
        k = make_k(data, frac, [1, 0, 0])
        d0_k = build_d0_bloch(V, E, k, L_vec, shifts)
        d1_k = build_d1_bloch_exact(V, E, F, k, L_vec, d0_k)
        norm = np.linalg.norm(d1_k @ d0_k)

        rank_d0 = np.linalg.matrix_rank(d0_k, tol=1e-10)
        K, M = build_K_M(d1_k, star1, star2)
        eigs = np.real(eigh(K, M, eigvals_only=True))
        n_zero = count_zeros(eigs)

        print(f"         frac={frac:.3f}: ‖d₁d₀‖={norm:.2e}, "
              f"rank(d₀)={rank_d0}, n_zero={n_zero}")
        # Exactness should hold everywhere
        assert norm < 1e-12, f"Exactness fails at frac={frac}"

    # Note: at k=k_BZ (frac=1.0), rank(d₀) drops by 1 (one Bloch phase becomes
    # trivial: e^{ik_j L_j}=1 along that axis). n_zero rises above V because
    # ker(d₁) gains additional harmonic modes (H¹ ≠ 0 at TRIM points).
    # This is expected from Künneth: at TRIM, some factors of S¹ have trivial
    # local coefficients → nontrivial cohomology. NOT a bug.
    print("         Note: n_zero > V at k_BZ is expected (TRIM cohomology)")
    print("         Exactness maintained at BZ boundary  [PASS]")


def main():
    print("=" * 60)
    print("GROUP 3: The fix — recurrence construction (§4)")
    print("=" * 60)

    structs = build_structures()
    test_exactness_all_structures(structs)
    test_holonomy(structs)
    test_canonical_K_permutations(structs)
    test_gauge_invariance(structs)
    test_kernel_equals_V(structs)
    test_rank_chain(structs)
    test_perturbation_stability(structs)
    test_d2_exactness(structs)
    test_cohomology(structs)
    test_cohomology_trim(structs)
    test_iff_completeness(structs)
    test_bz_boundary(structs)

    print("\n" + "-" * 60)
    print("GROUP 3: ALL PASS")
    print("-" * 60)


if __name__ == '__main__':
    main()
