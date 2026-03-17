"""
GROUP 1 — Setup validation (§2).

Verify test structures are correctly built and DEC operators are well-formed.

Tests:
  T1.1  d₀ has correct shape (E×V), entries in {-1, 0, +1}
  T1.2  d₁ topological: d₁d₀ = 0 at k=0
  T1.3  ⋆₁, ⋆₂ strictly positive (Voronoi structures only; SC has no dual)
  T1.4  Structure counts V, E, F match expected
  T1.5  Scalar Laplacian K₀ has no spurious zeros at k≠0 (3 directions)

Note: SC N=3 has no Voronoi dual → T1.3, T1.5 skip SC (documented, not silent).
      SC is tested for topology (T1.1, T1.2, T1.4) and for Bloch operators in Groups 2-5.
      Random Voronoi (3 seeds) included in build_structures().

Expected output (verified 17 Mar 2026, macOS, Python 3, NumPy/SciPy):

  ============================================================
  GROUP 1: Setup validation (§2)
  ============================================================

    T1.4  Structure counts
           Kelvin N=2: V=96 E=192 F=112  [PASS]
           C15 N=1: V=136 E=272 F=160  [PASS]
           WP N=1: V=46 E=92 F=54  [PASS]
           SC N=3: V=27 E=81 F=81  [PASS]
           Voronoi seed=42: V=331 E=662 F=381  [built]
           Voronoi seed=137: V=347 E=694 F=397  [built]
           Voronoi seed=999: V=354 E=708 F=404  [built]

    T1.1  d₀ shape and entries
           Kelvin N=2: (192×96), entries {-1,0,+1}  [PASS]
           C15 N=1: (272×136), entries {-1,0,+1}  [PASS]
           WP N=1: (92×46), entries {-1,0,+1}  [PASS]
           SC N=3: (81×27), entries {-1,0,+1}  [PASS]
           Voronoi seed=42: (662×331), entries {-1,0,+1}  [PASS]
           Voronoi seed=137: (694×347), entries {-1,0,+1}  [PASS]
           Voronoi seed=999: (708×354), entries {-1,0,+1}  [PASS]

    T1.2  d₁d₀ = 0 at k=0 (topological)
           Kelvin N=2: ||d₁d₀|| = 0.00e+00  [PASS]
           C15 N=1: ||d₁d₀|| = 0.00e+00  [PASS]
           WP N=1: ||d₁d₀|| = 0.00e+00  [PASS]
           SC N=3: ||d₁d₀|| = 0.00e+00  [PASS]
           Voronoi seed=42: ||d₁d₀|| = 0.00e+00  [PASS]
           Voronoi seed=137: ||d₁d₀|| = 0.00e+00  [PASS]
           Voronoi seed=999: ||d₁d₀|| = 0.00e+00  [PASS]

    T1.3  Hodge stars positive (Voronoi only; SC has no dual)
           Kelvin N=2: ⋆₁ min=4.0000, ⋆₂ min=0.6667  [PASS]
           C15 N=1: ⋆₁ min=1.1111, ⋆₂ min=1.4876  [PASS]
           WP N=1: ⋆₁ min=1.6000, ⋆₂ min=0.8205  [PASS]
           SC N=3: no Voronoi dual — tested for topology only (T1.1, T1.2, T1.4)
           Voronoi seed=42: ⋆₁ min=0.0630, ⋆₂ min=0.1636  [PASS]
           Voronoi seed=137: ⋆₁ min=0.0243, ⋆₂ min=0.1236  [PASS]
           Voronoi seed=999: ⋆₁ min=0.0452, ⋆₂ min=0.1042  [PASS]

    T1.5  Scalar Laplacian: no spurious zeros at k≠0 (3 directions)
           Kelvin N=2: [100] min=3.29e-02, [110] min=3.29e-02, [111] min=3.29e-02  [PASS]
           C15 N=1: [100] min=1.16e-02, [110] min=1.16e-02, [111] min=1.16e-02  [PASS]
           WP N=1: [100] min=3.43e-02, [110] min=3.43e-02, [111] min=3.43e-02  [PASS]
           SC N=3: no Voronoi dual — tested for topology only (T1.1, T1.2, T1.4)
           Voronoi seed=42: [100] min=4.76e-03, [110] min=4.76e-03, [111] min=4.76e-03  [PASS]
           Voronoi seed=137: [100] min=4.54e-03, [110] min=4.54e-03, [111] min=4.54e-03  [PASS]
           Voronoi seed=999: [100] min=4.46e-03, [110] min=4.45e-03, [111] min=4.45e-03  [PASS]

  ------------------------------------------------------------
  GROUP 1: ALL PASS (1.3s)
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
from physics.gauge_bloch import compute_edge_shifts, build_d0_bloch
from core_math.builders.solids_periodic import build_sc_supercell_periodic
from core_math.operators.incidence import build_d0, build_d1


# ── Expected counts ──────────────────────────────────────────────────────

EXPECTED = {
    'Kelvin N=2':  {'V': 96,  'E': 192, 'F': 112},
    'C15 N=1':     {'V': 136, 'E': 272, 'F': 160},
    'WP N=1':      {'V': 46,  'E': 92,  'F': 54},
    'SC N=3':      {'V': 27,  'E': 81,  'F': 81},
}

# Structures with Voronoi dual (have cell_centers → support Hodge stars)
HAS_VORONOI_DUAL = {'Kelvin N=2', 'C15 N=1', 'WP N=1'}


def wrap_sc(N):
    """Wrap SC tuple into dict. SC has no Voronoi dual."""
    verts, edges, faces, cells = build_sc_supercell_periodic(N=N)
    a = 2.0
    L = a * N
    return {
        'V': verts, 'E': edges, 'F': faces, 'C': cells,
        'L': L, 'L_vec': np.array([L, L, L]),
    }


def build_random_voronoi(seed, n_cells=50, L=4.0):
    """Build random periodic Voronoi. Returns (name, data) or None on failure."""
    np.random.seed(seed)
    points = np.random.uniform(0, L, size=(n_cells, 3))
    try:
        data = build_foam_with_dual_info(points, L)
        return f'Voronoi seed={seed}', data
    except Exception as e:
        print(f"  WARNING: Voronoi seed={seed} build failed: {e}")
        return None


def build_structures():
    """Build all test structures including random Voronoi."""
    structs = {}
    structs['Kelvin N=2'] = build_kelvin_with_dual_info(N=2)
    structs['C15 N=1'] = build_c15_with_dual_info(N=1)
    structs['WP N=1'] = build_wp_with_dual_info(N=1)
    structs['SC N=3'] = wrap_sc(N=3)

    for seed in [42, 137, 999]:
        result = build_random_voronoi(seed)
        if result:
            name, data = result
            structs[name] = data
            HAS_VORONOI_DUAL.add(name)
    return structs


def has_dual(name):
    return name in HAS_VORONOI_DUAL


# ── Tests ────────────────────────────────────────────────────────────────

def test_structure_counts(structs):
    """T1.4: V, E, F counts match expected values."""
    print("\n  T1.4  Structure counts")
    for name, data in structs.items():
        nV, nE, nF = len(data['V']), len(data['E']), len(data['F'])
        if name in EXPECTED:
            exp = EXPECTED[name]
            ok = (nV == exp['V'] and nE == exp['E'] and nF == exp['F'])
            status = 'PASS' if ok else 'FAIL'
            print(f"         {name}: V={nV} E={nE} F={nF}  [{status}]")
            assert ok, f"{name}: counts mismatch"
        else:
            # Random Voronoi — just report, no expected values
            print(f"         {name}: V={nV} E={nE} F={nF}  [built]")


def test_euler_characteristic(structs):
    """T1.6: Euler characteristic χ = V - E + F - C = 0 on T³."""
    print("\n  T1.6  Euler characteristic χ = 0 on T³")
    for name, data in structs.items():
        nV, nE, nF = len(data['V']), len(data['E']), len(data['F'])
        if 'C' in data:
            nC = len(data['C'])
            chi = nV - nE + nF - nC
            assert chi == 0, f"{name}: χ = {chi} ≠ 0"
            print(f"         {name}: V={nV} E={nE} F={nF} C={nC}, χ={chi}  [PASS]")
        else:
            # No cells — check V - E + F (should be 0 for surface of T³ complex)
            chi_2 = nV - nE + nF
            print(f"         {name}: V={nV} E={nE} F={nF}, V-E+F={chi_2}  [no cells]")


def test_d0_shape_and_entries(structs):
    """T1.1: d₀ has shape (E×V) with entries in {-1, 0, +1}."""
    print("\n  T1.1  d₀ shape and entries")
    for name, data in structs.items():
        V, E = data['V'], data['E']
        nV, nE = len(V), len(E)
        d0 = build_d0(V, E)
        assert d0.shape == (nE, nV), f"{name}: d0 shape {d0.shape} != ({nE}, {nV})"
        vals = set(np.unique(d0))
        assert vals <= {-1, 0, 1}, f"{name}: d0 has values {vals}"
        for i in range(nE):
            row = d0[i, :]
            assert np.sum(row == 1) == 1 and np.sum(row == -1) == 1, \
                f"{name}: d0 row {i} malformed"
        print(f"         {name}: ({nE}×{nV}), entries {{-1,0,+1}}  [PASS]")


def test_d1d0_at_k0(structs):
    """T1.2: d₁d₀ = 0 at k=0 (topological, no Bloch phases)."""
    print("\n  T1.2  d₁d₀ = 0 at k=0 (topological)")
    for name, data in structs.items():
        V, E, F = data['V'], data['E'], data['F']
        d0 = build_d0(V, E)
        d1 = build_d1(V, E, F)
        norm = np.linalg.norm(d1 @ d0)
        assert norm < 1e-14, f"{name}: ||d1*d0|| = {norm:.2e} at k=0"
        print(f"         {name}: ||d₁d₀|| = {norm:.2e}  [PASS]")


def test_hodge_stars_positive(structs):
    """T1.3: ⋆₁ and ⋆₂ strictly positive (Voronoi structures only; SC has no dual)."""
    print("\n  T1.3  Hodge stars positive (Voronoi only; SC has no dual)")
    for name, data in structs.items():
        if not has_dual(name):
            print(f"         {name}: no Voronoi dual — tested for topology only (T1.1, T1.2, T1.4)")
            continue
        star1, star2 = build_hodge_stars_voronoi(data)
        assert np.all(star1 > 0), f"{name}: star1 has non-positive entries"
        assert np.all(star2 > 0), f"{name}: star2 has non-positive entries"
        print(f"         {name}: ⋆₁ min={np.min(star1):.4f}, ⋆₂ min={np.min(star2):.4f}  [PASS]")


def test_scalar_laplacian_no_spurious(structs):
    """T1.5: Scalar Laplacian K₀ has no zero eigenvalues at k≠0 (3 directions)."""
    print("\n  T1.5  Scalar Laplacian: no spurious zeros at k≠0 (3 directions)")
    directions = {
        '[100]': [1.0, 0.0, 0.0],
        '[110]': [1.0, 1.0, 0.0],
        '[111]': [1.0, 1.0, 1.0],
    }
    for name, data in structs.items():
        if not has_dual(name):
            print(f"         {name}: no Voronoi dual — tested for topology only (T1.1, T1.2, T1.4)")
            continue
        V, E = data['V'], data['E']
        L_vec = data['L_vec']
        star1, _ = build_hodge_stars_voronoi(data)
        shifts = compute_edge_shifts(V, E, L_vec)

        mins = {}
        for d_label, d_vec in directions.items():
            d_norm = np.array(d_vec) / np.linalg.norm(d_vec)
            k = 2 * np.pi / data['L'] * 0.1 * d_norm
            d0_k = build_d0_bloch(V, E, k, L_vec, shifts)
            K0 = d0_k.conj().T @ np.diag(star1) @ d0_k
            K0 = 0.5 * (K0 + K0.conj().T)
            eigs = np.real(eigh(K0, eigvals_only=True))
            n_zero = np.sum(np.abs(eigs) < 1e-10)
            assert n_zero == 0, f"{name} {d_label}: K₀ has {n_zero} zeros"
            mins[d_label] = np.min(eigs)

        parts = ", ".join(f"{d} min={mins[d]:.2e}" for d in directions)
        print(f"         {name}: {parts}  [PASS]")


def main():
    print("=" * 60)
    print("GROUP 1: Setup validation (§2)")
    print("=" * 60)

    structs = build_structures()
    test_structure_counts(structs)
    test_euler_characteristic(structs)
    test_d0_shape_and_entries(structs)
    test_d1d0_at_k0(structs)
    test_hodge_stars_positive(structs)
    test_scalar_laplacian_no_spurious(structs)

    print("\n" + "-" * 60)
    print("GROUP 1: ALL PASS")
    print("-" * 60)


if __name__ == '__main__':
    main()
