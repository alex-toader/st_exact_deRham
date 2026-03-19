# Tests Map — st_exact_deRham (ETNA version)

**Date:** 19 Mar 2026
**Target:** Electronic Transactions on Numerical Analysis (ETNA)
**Total:** 7 files, 53 tests
**Runtime:** ~70s full suite
**Run:** `cd st_exact_deRham && PYTHONPATH=src python3 tests/test_0X_*.py`

---

## Mapping to ETNA paper sections

| Test file | Paper §  | What it verifies |
|-----------|----------|------------------|
| test_01_setup | §2 (Setup) | Cell complex, incidence, Hodge stars |
| test_02_standard_fails | §3 (Failure) | Prop 1, Cor 1 |
| test_03_recurrence | §4 (Construction) | Thm 1, Lemma, Prop 2 |
| test_04_spectral_gains | §5 (Numerics) | Table 1, Figure 1, convergence |
| test_05_inexactness | — (cut from ETNA) | Rank-pollution, hybridization, trace |
| test_06_voronoi | §6 (Voronoi remark) | G=Vol·I, c²=1 |
| test_07_curvature | — (cut from ETNA) | Curvature structure |

Tests 05 and 07 verify claims cut from the ETNA note. They remain in
the test suite for completeness and for future papers (A2).

---

## test_01_setup.py — §2 Setup (6 tests)

| Test | Paper claim | Key result |
|------|---------------|------------|
| T1.1 | d₀ shape, entries {-1,0,+1} | PASS all 7 structures |
| T1.2 | d₁d₀ = 0 at k=0 (topological) | ‖d₁d₀‖ = 0.00 |
| T1.3 | ⋆₁, ⋆₂ > 0 (Voronoi) | min > 0 all |
| T1.4 | V,E,F counts | PASS |
| T1.5 | No spurious zeros at k≠0 | PASS 3 dirs |
| T1.6 | Euler χ = 0 on T³ | PASS |

## test_02_standard_fails.py — §3 Failure (12 tests)

| Test | Paper claim | Key result |
|------|---------------|------------|
| T2.1 | ‖d₁_std d₀‖ > 1 | 5.07–11.68 → Prop 1 |
| T2.2 | Direction-dependent pollution | [100]<[110]<[111] → Prop 1 |
| T2.3 | n_spur > 0 | 3–19 → Prop 1 |
| T2.4 | Intra-face contradictions | 37/112 → Cor 1 |
| T2.5 | Per-face 1D null space | 112/112 → Cor 1 |
| T2.6 | ‖d₁_std d₀‖ ~ O(k) | slope 0.96 |
| T2.7 | ker(K_std) < V | deficit 3–19 |
| T2.8 | rank(d₁_std) direction-dependent | 102/108/110 vs 96 |
| T2.9 | k=0 control: exact = standard | PASS |
| T2.10 | n_spur scaling ~ N^2.3 | N^2.29–2.37 |
| T2.11 | Per-edge optimization min = 1.2564 | ratio 10¹⁵ → Cor 1 |
| T2.12 | Phase sensitivity O(ε) | ε-linear, isolated → Thm 1 remark |

## test_03_recurrence.py — §4 Construction (12 tests)

| Test | Paper claim | Key result |
|------|---------------|------------|
| T3.1 | ‖d₁_exact d₀‖ < 10⁻¹⁴ | max 9.53e-16 → **Thm 1** |
| T3.2 | Holonomy H_f = 1 | max |H-1| = 2.3e-16 → **Lemma** |
| T3.3 | K canonical (permutations) | ‖ΔK‖ < 5.2e-15 → Thm 1 |
| T3.4 | K gauge-invariant | max |Δeig| < 10⁻¹² → Thm 1 |
| T3.5 | n_zero = V all k | 32 k × 4 structs → **Prop 2** |
| T3.6 | rank chain im(d₀)=ker(d₁)=ker(K) | grad < 1.3e-13 → Prop 2 |
| T3.7 | Exactness under vertex perturbation | ε=0..50% → Thm 1 |
| T3.8 | d₂d₁ = 0 (full complex) | 10⁻¹⁵ → §4 last paragraph |
| T3.9 | β(Γ)=(1,3,3,1), β(k≠0)=(0,0,0,0) | → §4 last paragraph |
| T3.10 | BZ boundary: exactness OK | TRIM points → Prop 2 |
| T3.11 | nullity = nF (iff completeness) | → Thm 1 uniqueness |
| T3.12 | Künneth at TRIM | β₁=0 at X,M,R → Prop 2 |

## test_04_spectral_gains.py — §5 Numerics (7 tests)

| Test | Paper claim | Key result |
|------|---------------|------------|
| T4.1 | Hodge splitting | exact 0%, std 73–92% → §5 |
| T4.2 | Gradient overlap exact | max 1.34e-13 → §5 |
| T4.3 | Spurious modes gradient | mean 0.88–0.94 → §5 |
| T4.4 | c² convergence Kelvin | p=2.00 → §5 |
| T4.5 | c² convergence SC | p=2.00 → §5 |
| T4.7 | Universality 4 structures | PASS → §5 |
| T4.8 | Random Voronoi 10 seeds | 10/10 → §5 |

## test_05_inexactness.py — cut from ETNA (3 tests)

| Test | Claim (cut from ETNA) | Key result |
|------|----------------|------------|
| T5.1 | rank excess = n_spur | 6=6, 9=9, 3=3 |
| T5.2 | hybridization → random baseline | 0.514, 0.524, 0.513 |
| T5.3 | tr(M⁻¹K) conserved | rel_diff < 10⁻¹⁶ |

## test_06_voronoi.py — §6 Voronoi remark (8 tests)

| Test | Paper claim | Key result |
|------|---------------|------------|
| T6.1 | G=Vol·I cubic | ‖G/Vol−I‖ < 10⁻¹⁵ → Remark §6 |
| T6.2 | G=Vol·I random | < 4.3e-12 → Remark §6 |
| T6.3 | H=Vol·I cubic | < 7.8e-16 → Remark §6 |
| T6.4 | c²=1 cubic (exact) | 0.9993–0.9997 → Remark §6 |
| T6.5 | c²=1 random Voronoi | 0.9995–0.9996 → Remark §6 |
| T6.6 | c²_std ≠ 1 | 0.52–1.55 → Remark §6 |
| T6.7 | perturbed stars → c≠1 | Δc²=0.012/0.024 |
| T6.8 | G=Vol·I near-degenerate | err=3.5e-12 |

## test_07_curvature.py — cut from ETNA (5 tests)

| Test | Claim (cut from ETNA) | Key result |
|------|-------------------|------------|
| T7.1 | ratio unimodular | max dev 2.2e-16 |
| T7.2 | n_curv topological | k-independent |
| T7.3 | flux ∝ k | linear, CV=0 |
| T7.4 | n_curv > n_spur | cancellation |
| T7.5 | flux structure-dependent | integers |
