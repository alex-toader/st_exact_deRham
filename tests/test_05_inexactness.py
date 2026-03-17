"""
GROUP 5 — Consequences of inexactness (§6, condensed).

3 results explaining WHY standard DEC fails spectrally.
These are consequences of Prop 1, explaining the numerics in §5.

Tests:
  T5.1  Rank pollution = rank excess (algebraic identity)
  T5.2  Hybridization: spurious overlap → nV/nE ≈ 0.5 (random baseline)
  T5.3  Trace conservation: Σ(eigenvalue shifts) + Σ(spurious eigs) = 0

Source: W3 files 3, 5, 6 from st_bloch_exactness.

Expected output (verified 17 Mar 2026):

  ============================================================
  GROUP 5: Consequences of inexactness (§6)
  ============================================================

    T5.1  Rank pollution = rank excess (algebraic identity)
           Kelvin N=2: rank_excess = rank_std(102) − rank_exact(96) = 6 = n_spur  [PASS]
           C15 N=1: rank_excess = rank_std(145) − rank_exact(136) = 9 = n_spur  [PASS]
           WP N=1: rank_excess = rank_std(49) − rank_exact(46) = 3 = n_spur  [PASS]
           Kelvin direction scan:
             [100]: Δrank=6 = n_spur=6  [PASS]
             [110]: Δrank=12 = n_spur=12  [PASS]
             [111]: Δrank=14 = n_spur=14  [PASS]

    T5.2  Hybridization at random baseline nV/nE
           Kelvin N=2: exact(g=1.000, p=0.000), spur=0.514, baseline=0.500  [PASS]
           C15 N=1: exact(g=1.000, p=0.000), spur=0.524, baseline=0.500  [PASS]
           WP N=1: exact(g=1.000, p=0.000), spur=0.513, baseline=0.500  [PASS]
           Universal: spurious at random baseline on all structures  [PASS]

    T5.3  Trace conservation: tr(M⁻¹K) exact = standard
           Kelvin N=2: tr(M⁻¹K)_ex=160.0000, tr(M⁻¹K)_std=160.0000, rel_diff=1.78e-16  [PASS]
           C15 N=1: tr(M⁻¹K)_ex=1047.7792, tr(M⁻¹K)_std=1047.7792, rel_diff=0.00e+00  [PASS]
           WP N=1: tr(M⁻¹K)_ex=174.3732, tr(M⁻¹K)_std=174.3732, rel_diff=0.00e+00  [PASS]
           Multi-k verification (Kelvin, 5 k-points):
             max rel_diff across 5 k-points: 1.78e-16  [PASS]

  ------------------------------------------------------------
  GROUP 5: ALL PASS (1.7s)
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
from physics.bloch import build_d1_bloch_standard, compute_edge_crossings, build_edge_lookup


# ── Helpers ──────────────────────────────────────────────────────────────

def build_structures():
    return {
        'Kelvin N=2': build_kelvin_with_dual_info(N=2),
        'C15 N=1':    build_c15_with_dual_info(N=1),
        'WP N=1':     build_wp_with_dual_info(N=1),
    }


def make_k(data, frac, direction):
    d = np.array(direction, dtype=float)
    d /= np.linalg.norm(d)
    return 2 * np.pi / data['L'] * frac * d


def build_both(data, k):
    """Build exact and standard d₁, return all operators."""
    V, E, F = data['V'], data['E'], data['F']
    L, L_vec = data['L'], data['L_vec']
    shifts = compute_edge_shifts(V, E, L_vec)
    crossings = compute_edge_crossings(V, E, L)
    edge_lookup = build_edge_lookup(E, crossings)
    star1, star2 = build_hodge_stars_voronoi(data)

    d0_k = build_d0_bloch(V, E, k, L_vec, shifts)
    d1_ex = build_d1_bloch_exact(V, E, F, k, L_vec, d0_k)
    d1_std = build_d1_bloch_standard(V, E, F, L, k, edge_lookup, crossings)

    M = np.diag(star1)

    def build_K(d1):
        K = d1.conj().T @ np.diag(star2) @ d1
        return 0.5 * (K + K.conj().T)

    return {
        'd0_k': d0_k, 'd1_ex': d1_ex, 'd1_std': d1_std,
        'K_ex': build_K(d1_ex), 'K_std': build_K(d1_std),
        'M': M, 'star1': star1, 'star2': star2,
        'nV': len(V), 'nE': len(E),
    }


def grad_overlap(u, d0, M):
    """Gradient fraction: ||P_grad u||²_M / ||u||²_M."""
    Md0 = M @ d0
    G_inv = np.linalg.inv(d0.conj().T @ Md0)
    c = G_inv @ (d0.conj().T @ M @ u)
    Pu = d0 @ c
    return np.real(Pu.conj() @ M @ Pu) / np.real(u.conj() @ M @ u)


# ── Tests ────────────────────────────────────────────────────────────────

def test_rank_pollution_identity(structs):
    """T5.1: rank(d₁_std) − rank(d₁_exact) = n_spur. Algebraic identity.

    This is not empirical — it follows from both operators having the same
    sparsity pattern (same nonzero positions) but different phases.
    The rank difference directly gives the pollution count.
    """
    print("\n  T5.1  Rank pollution = rank excess (algebraic identity)")
    for name, data in structs.items():
        k = make_k(data, 0.10, [1, 0, 0])
        ops = build_both(data, k)

        rank_ex = np.linalg.matrix_rank(ops['d1_ex'], tol=1e-10)
        rank_std = np.linalg.matrix_rank(ops['d1_std'], tol=1e-10)
        rank_excess = rank_std - rank_ex

        # Count spurious from eigenvalues
        eigs_std = np.real(eigh(ops['K_std'], ops['M'], eigvals_only=True))
        thresh = max(np.max(np.abs(eigs_std)) * 1e-12, 1e-10)
        n_zero_std = int(np.sum(np.abs(eigs_std) < thresh))
        n_spur = ops['nV'] - n_zero_std

        assert rank_excess == n_spur, \
            f"{name}: rank excess {rank_excess} ≠ n_spur {n_spur}"
        print(f"         {name}: rank_excess = rank_std({rank_std}) − "
              f"rank_exact({rank_ex}) = {rank_excess} = n_spur  [PASS]")

    # Verify on multiple directions to show it's not accidental
    data = structs['Kelvin N=2']
    print(f"         Kelvin direction scan:")
    for d_label, d_vec in [('[100]', [1,0,0]), ('[110]', [1,1,0]), ('[111]', [1,1,1])]:
        k = make_k(data, 0.10, d_vec)
        ops = build_both(data, k)
        rank_ex = np.linalg.matrix_rank(ops['d1_ex'], tol=1e-10)
        rank_std = np.linalg.matrix_rank(ops['d1_std'], tol=1e-10)
        eigs_std = np.real(eigh(ops['K_std'], ops['M'], eigvals_only=True))
        thresh = max(np.max(np.abs(eigs_std)) * 1e-12, 1e-10)
        n_spur = ops['nV'] - int(np.sum(np.abs(eigs_std) < thresh))
        assert rank_std - rank_ex == n_spur
        print(f"           {d_label}: Δrank={rank_std - rank_ex} = n_spur={n_spur}  [PASS]")


def test_hybridization_random_baseline(structs):
    """T5.2: Spurious mode overlap → nV/nE (random baseline). Universal.

    On exact DEC: gauge modes have overlap 1.0, physical modes 0.0.
    On standard DEC: spurious modes have overlap ≈ nV/nE ≈ 0.5 — the
    value expected for a RANDOM vector in the edge space. This means
    standard DEC spurious modes have literally no preferred direction
    in the gradient/curl decomposition.
    """
    print("\n  T5.2  Hybridization at random baseline nV/nE")
    # k chosen generic (non-axis, non-diagonal) to avoid symmetry-protected directions
    # that could mask hybridization.
    k_frac = np.array([0.17, 0.23, 0.31])

    for name, data in structs.items():
        V, E, F = data['V'], data['E'], data['F']
        L, L_vec = data['L'], data['L_vec']
        nV, nE = len(V), len(E)
        star1, star2 = build_hodge_stars_voronoi(data)
        shifts = compute_edge_shifts(V, E, L_vec)
        crossings = compute_edge_crossings(V, E, L)
        edge_lookup = build_edge_lookup(E, crossings)
        M = np.diag(star1)

        k = 2 * np.pi / L * k_frac
        d0_k = build_d0_bloch(V, E, k, L_vec, shifts)
        d1_ex = build_d1_bloch_exact(V, E, F, k, L_vec, d0_k)
        d1_std = build_d1_bloch_standard(V, E, F, L, k, edge_lookup, crossings)

        # Exact: eigensolve
        K_ex = d1_ex.conj().T @ np.diag(star2) @ d1_ex
        K_ex = 0.5 * (K_ex + K_ex.conj().T)
        evals_ex, evecs_ex = eigh(K_ex, M)

        # Standard: eigensolve
        K_std = d1_std.conj().T @ np.diag(star2) @ d1_std
        K_std = 0.5 * (K_std + K_std.conj().T)
        evals_std, evecs_std = eigh(K_std, M)

        # Classify standard modes
        n_gauge_std = int(np.sum(np.abs(evals_std) < 1e-8))
        n_spur = nV - n_gauge_std

        # Gradient overlap using EXACT d₀ as reference
        ex_gauge_ov = np.mean([grad_overlap(evecs_ex[:, i], d0_k, M) for i in range(nV)])
        ex_phys_ov = np.mean([grad_overlap(evecs_ex[:, i], d0_k, M) for i in range(nV, nE)])

        spur_ov = np.mean([grad_overlap(evecs_std[:, i], d0_k, M)
                           for i in range(n_gauge_std, nV)]) if n_spur > 0 else 0.0

        baseline = nV / nE

        assert abs(ex_gauge_ov - 1.0) < 1e-4, f"{name}: exact gauge overlap ≠ 1"
        assert abs(ex_phys_ov) < 1e-4, f"{name}: exact phys overlap ≠ 0"
        assert abs(spur_ov - baseline) < 0.04, \
            f"{name}: spurious overlap {spur_ov:.3f} not near baseline {baseline:.3f}"

        print(f"         {name}: exact(g={ex_gauge_ov:.3f}, p={ex_phys_ov:.3f}), "
              f"spur={spur_ov:.3f}, baseline={baseline:.3f}  [PASS]")

    print(f"         Universal: spurious at random baseline on all structures  [PASS]")


def test_trace_conservation(structs):
    """T5.3: tr(M⁻¹K) conserved between exact and standard.

    Total trace tr(M⁻¹K) = Σ_f ⋆₂[f] × Σ_{e∈∂f} ⋆₁[e]⁻¹ |d₁[f,e]|².
    Since |d₁[f,e]| = 1 for both constructions (same support, unimodular),
    the trace is identical. This is a topological invariant.
    Consequence: spurious eigenvalues must be compensated by shifts
    in physical eigenvalues. The total is conserved.
    """
    print("\n  T5.3  Trace conservation: tr(M⁻¹K) exact = standard")
    for name, data in structs.items():
        k = make_k(data, 0.10, [1, 0, 0])
        ops = build_both(data, k)

        # Method 1: matrix trace tr(M⁻¹K)
        M_inv = np.diag(1.0 / ops['star1'])
        tr_mat_ex = np.real(np.trace(M_inv @ ops['K_ex']))
        tr_mat_std = np.real(np.trace(M_inv @ ops['K_std']))

        # Method 2: Frobenius norm of d₁ (direct: |d₁[f,e]|=1 on same support)
        fr_ex = np.linalg.norm(ops['d1_ex'], 'fro')
        fr_std = np.linalg.norm(ops['d1_std'], 'fro')

        rel_diff_mat = abs(tr_mat_ex - tr_mat_std) / abs(tr_mat_ex)
        fr_diff = abs(fr_ex - fr_std)

        assert rel_diff_mat < 1e-12, \
            f"{name}: tr(M⁻¹K) differs by {rel_diff_mat:.2e}"
        assert fr_diff < 1e-12, \
            f"{name}: ||d₁||_F differs by {fr_diff:.2e}"

        print(f"         {name}: tr(M⁻¹K)_ex={tr_mat_ex:.4f}, tr(M⁻¹K)_std={tr_mat_std:.4f}, "
              f"rel_diff={rel_diff_mat:.2e}, ‖d₁‖_F diff={fr_diff:.2e}  [PASS]")

    # Multi-k verification: trace conservation is k-independent (topological)
    print(f"         Multi-k verification (Kelvin, 5 k-points):")
    data = structs['Kelvin N=2']
    k_dirs = [[1,0,0], [0,1,0], [1,1,0], [1,1,1], [0.17,0.23,0.31]]
    max_rel = 0
    for d in k_dirs:
        k = make_k(data, 0.10, d)
        ops = build_both(data, k)
        M_inv = np.diag(1.0 / ops['star1'])
        tr_ex = np.real(np.trace(M_inv @ ops['K_ex']))
        tr_std = np.real(np.trace(M_inv @ ops['K_std']))
        rel = abs(tr_ex - tr_std) / abs(tr_ex)
        max_rel = max(max_rel, rel)
    assert max_rel < 1e-12, f"Trace conservation fails at some k: {max_rel:.2e}"
    print(f"           max rel_diff across 5 k-points: {max_rel:.2e}  [PASS]")


def main():
    print("=" * 60)
    print("GROUP 5: Consequences of inexactness (§6)")
    print("=" * 60)

    structs = build_structures()
    test_rank_pollution_identity(structs)
    test_hybridization_random_baseline(structs)
    test_trace_conservation(structs)

    print("\n" + "-" * 60)
    print("GROUP 5: ALL PASS")
    print("-" * 60)


if __name__ == '__main__':
    main()
