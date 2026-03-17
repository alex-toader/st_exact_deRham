# st_exact_deRham 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19074708.svg)](https://doi.org/10.5281/zenodo.19074708)

Exactness-preserving discrete de Rham complexes under Bloch-periodic boundary conditions.

**Author:** Alexandru Toader (toader_alexandru@yahoo.com)

## Result

Standard Bloch-periodic DEC on unstructured polyhedral meshes breaks the discrete de Rham exact sequence: d₁(k)d₀(k) ≠ 0 for generic wave vector k. No per-edge phase assignment can restore exactness (the obstruction requires face-dependent phases, not edge-dependent). A face-boundary recurrence constructs the unique (up to per-face gauge) flat extension of the Bloch-twisted gradient to a cochain complex with d₁(k)d₀(k) = 0 for all k. The curl–curl operator K(k) is canonical and has gauge kernel dimension |V|.

Structural consequences of inexactness: the pollution count equals the rank excess of d₁ (algebraic identity), the Hodge decomposition degrades to random-baseline hybridization, and the eigenvalue trace is conserved (topological invariant). The curvature of the standard construction relative to the exact one has a rigid topological structure: integer flux coefficients, direction-dependent support, and partial cancellation (n_curv > n_spur).

On periodic Voronoi tessellations, exactness combined with metric isotropy G = H = Vol·I yields ω² = |k|² + O(|k|⁴).

Verified on SC cubic, Kelvin (BCC), C15 (Laves), Weaire-Phelan, and random Voronoi complexes.

## Paper

Submitted to SIAM Journal on Numerical Analysis.

- `paper/main.tex` — manuscript
- `paper/supplementary.tex` — supplementary material (Voronoi proofs)
- `paper/main.pdf` — compiled PDF
- `paper/supplementary.pdf` — compiled SM

## Tests

57 tests across 7 files (~70s). Each test file maps to a paper section.

See `tests/tests_map.md` for the complete inventory with per-test descriptions.

```
tests/
├── test_01_setup.py           (6 tests)   Structure validation, Euler χ, Hodge stars
├── test_02_standard_fails.py  (12 tests)  Exactness failure, Cor 1, scaling, sensitivity
├── test_03_recurrence.py      (14 tests)  Theorem 1, holonomy, uniqueness, kernel, cohomology, iff
├── test_04_spectral_gains.py  (8 tests)   Pollution, Hodge splitting, convergence, universality
├── test_05_inexactness.py     (3 tests)   Rank identity, hybridization, trace conservation
├── test_06_voronoi.py         (8 tests)   Metric isotropy, c²=1, necessity, near-degeneracy
├── test_07_curvature.py       (5 tests)   Curvature structure, flux, cancellation
└── tests_map.md                           Complete test inventory
```

## Source code

```
src/
├── physics/
│   ├── gauge_bloch.py        Exactness-preserving d₁(k) construction
│   ├── bloch_complex.py      d₂(k) construction (full cochain complex)
│   ├── hodge.py              Voronoi complex builders + Hodge stars
│   ├── bloch.py              Standard Bloch operators (for comparison)
│   └── constants.py          Physical constants
└── core_math/
    ├── operators/incidence.py    Topological d₀, d₁ (k=0)
    ├── builders/                 Periodic structure builders (SC, Kelvin, C15, WP)
    └── spec/                     Mesh contract + constants
```

The core contribution is in `src/physics/gauge_bloch.py`: the function `build_d1_bloch_exact()` implements the recurrence construction. Compare with `build_d1_bloch_standard()` in `src/physics/bloch.py` which breaks exactness at k ≠ 0.

## Running tests

```bash
python3 run_tests.py                    # all tests (~70s)
python3 run_tests.py recurrence voronoi # specific groups
python3 run_tests.py setup standard_fails recurrence spectral_gains inexactness voronoi curvature  # full suite explicit
```

## Requirements

- Python 3.9+
- NumPy, SciPy
- Matplotlib (for figure generation only)

## License

MIT
