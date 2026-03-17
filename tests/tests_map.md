# Tests Map — st_exact_deRham (Paper A1)

**Date:** 17 Mar 2026
**Total:** 6 files, 42 tests (38 Must + 4 Nice + 4 NEG + 1 CONTROL + 2 PARTIAL)
**Runtime:** ~35s full suite
**Gap analysis:** See `internal/gap_analysis.md` — no critical gaps.

---

## test_01_setup.py — Group 1: Setup validation (§2) — 6 tests, 1.3s

| Test | What | Structures | Key result |
|------|------|-----------|------------|
| T1.1 | d₀ shape (E×V), entries {-1,0,+1} | All 7 | PASS |
| T1.2 | d₁d₀ = 0 at k=0 (topological) | All 7 | ‖d₁d₀‖ = 0.00 |
| T1.3 | ⋆₁, ⋆₂ > 0 (Voronoi only) | 6/7 (SC skip) | ⋆₁ min=0.024, ⋆₂ min=0.10 |
| T1.4 | V,E,F counts match expected | 4 cubic + 3 Voronoi | PASS |
| T1.5 | K₀ no spurious zeros (3 dirs) | 6/7 (SC skip) | n_zero = 0 all dirs |
| T1.6 | Euler χ = 0 on T³ | SC: χ=0, foams: V-E+F=C | PASS |

**Structures:** Kelvin N=2, C15 N=1, WP N=1, SC N=3, Voronoi seeds 42/137/999.

---

## test_02_standard_fails.py — Group 2: The problem (§3) — 12 tests, 33s

| Test | What | Key numbers | Claim |
|------|------|------------|-------|
| T2.1 | ‖d₁_std d₀‖ > 1 all structures | 5.07–11.68 (incl. random Voronoi) | Prop 1 |
| T2.2 | Pollution direction-dependent (3 structs) | [100]<[110]<[111] universal | Prop 1 |
| T2.3 | n_spur > 0 all structures | 3–19 | Prop 1 |
| T2.4 | Intra-face contradictions | 37/112 (Kelvin) | Cor 1 |
| T2.5 | Per-face 1D null space | 112/112 | Cor 1 |
| T2.6 | ‖d₁_std d₀‖ ~ O(k) | slope 0.96 (k≤0.10) | — |
| T2.7 | **NEG:** ker(K_std) < V | deficit 3–19 | Prop 3 contrast |
| T2.8 | **NEG:** rank(d₁_std) ≠ exact AND direction-dependent | 102/108/110 vs 96 | cohomology |
| T2.9 | **CONTROL:** k=0 exact = standard | |d₁_ex(0)| = |d₁_top| | Prop 1 scope |
| T2.10 | n_spur scaling ~ N^2.3 (extensive, surface) | [100]: N^2.29, [111]: N^2.37 | §3.3 |
| T2.11 | Per-edge optimization min = 1.2564 (ratio 10¹⁵) | unique global min, 10 starts | Cor 1 numeric |
| T2.12 | Phase sensitivity O(ε), exact is isolated | ‖d₁d₀‖/ε = 33.7, CV=0.001 | §4 remark |

---

## test_03_recurrence.py — Group 3: The fix (§4) — 12 tests, 1.5s

| Test | What | Key numbers | Claim |
|------|------|------------|-------|
| T3.1 | ‖d₁_exact d₀‖ < 10⁻¹⁴ all k | max 9.53e-16 (32 k-points) | **Thm 1** |
| T3.2 | Holonomy H_f = 1 | max |H-1| = 2.3e-16 (all faces) | Lem 1 |
| T3.3 | K canonical under permutations | ‖ΔK‖ < 5.2e-15 | Prop 2 |
| T3.4 | K gauge-invariant | max |Δeig| < 10⁻¹² (5 gauges) | Prop 2 |
| T3.5 | n_zero = V all k | 32 k-points × 4 structs | Prop 3 |
| T3.6 | rank chain im(d₀)=ker(d₁)=ker(K) | max grad = 1.3e-13 | Prop 3 |
| T3.7 | Exactness under ε=0..50% | ‖d₁d₀‖ = 7.86e-16 at all ε | Thm 1 |
| T3.8 | d₂(k)d₁(k) = 0 (full complex) | 10⁻¹⁵ all 4 structures | §4.6 |
| T3.9 | β(Γ)=(1,3,3,1), β(k≠0)=(0,0,0,0) | all 4 structures | Remark |
| T3.10 | BZ boundary: exactness OK | n_zero=98 at k_BZ (TRIM) | Nice |
| T3.11 | iff: nullity = nF (solution unique up to 1D/face) | Kelvin 112, C15 160, WP 54, SC 81 | Master thm (ii) |
| T3.12 | Künneth at TRIM: β₁=0 at X,M,R; β₁=3 at BZ_x≡Γ | all 4 TRIM match | §4.6 |

---

## test_04_spectral_gains.py — Group 4: What you gain (§5) — 8 tests, 27.9s

| Test | What | Key numbers | Result |
|------|------|------------|--------|
| T4.1 | Hodge splitting | exact 0%, standard 73–92% mixed | R1 |
| T4.2 | Gradient overlap exact | max 1.34e-13 | R2 |
| T4.3 | Spurious modes gradient | mean 0.879–0.937 | R2 |
| T4.4 | c² convergence Kelvin + std | exact p=2.00, std mean|c²-1|=0.75 | R3 |
| T4.5 | c² convergence SC | p=2.00 (convergence rate, not isotropy) | R3 |
| T4.7 | Universality 4 structures | all pass | — |
| T4.8 | Random Voronoi 10 seeds | 10/10 ‖d₁d₀‖ < 10⁻¹² | — |

---

## test_05_inexactness.py — Group 5: Inexactness consequences (§6) — 3 tests, 1.7s

| Test | What | Key numbers | Result |
|------|------|------------|--------|
| T5.1 | rank excess = n_spur (identity) | 6=6, 9=9, 3=3 + 3 dirs | R5 |
| T5.2 | Spur overlap → nV/nE | 0.514, 0.524, 0.513 (baseline 0.500) | R6 |
| T5.3 | tr(M⁻¹K) conserved + multi-k | rel_diff < 1.78e-16, 5 k-dirs | R7 |

---

## test_06_voronoi.py — Group 6: Voronoi optimality (§7) — 7 tests, 3.1s

| Test | What | Key numbers | Claim |
|------|------|------------|-------|
| T6.1 | G=Vol·I cubic | ‖G/Vol−I‖ < 8.9e-16 | Prop 4 |
| T6.2 | G=Vol·I random (5 seeds) | ‖G/Vol−I‖ < 4.3e-12 | Prop 4 |
| T6.3 | H=Vol·I cubic | ‖H/Vol−I‖ < 7.8e-16 | Prop 4 |
| T6.4 | c²=1 cubic (exact) | 0.9993–0.9997 | R8 |
| T6.5 | c²=1 random Voronoi | 0.9995–0.9996 | R8 |
| T6.6 | **NEG:** c²_std ≠ 1 (spurious) | 0.52–1.55 | R8 contrast |
| T6.7 | **NEG:** perturbed ⋆₁,⋆₂ → c≠1 | Δc²=0.012 (⋆₁), 0.024 (⋆₂) | — |
| T6.8 | G=Vol·I on near-degenerate Voronoi | err=3.5e-12, min⋆₁=3.5e-5 | Prop 4 robust |

---

## test_07_curvature.py — Group 7: Curvature structure (§4/§6) — 5 tests, 0.5s

| Test | What | Key numbers | Result |
|------|------|------------|--------|
| T7.1 | Ratio d₁_std/d₁_exact is unimodular | max dev 2.2e-16 | pure phase |
| T7.2 | n_curv topological (k-magnitude-independent) | 24/40/51 constant | topological |
| T7.3 | Flux ∝ k (linear, not quantized) | [100]: coeff=6.0 exact, CV=0 | linear |
| T7.4 | n_curv > n_spur (curvature partially cancels) | 24>6, 40>12, 51>14 | cancellation |
| T7.5 | Flux coefficient structure-dependent | Kelvin=6, C15=21, WP=5 (integers!) | structural |

---

## Summary

| Group | File | § | Tests | PASS | NEG | CONTROL | Time |
|-------|------|---|-------|------|-----|---------|------|
| 1 | test_01_setup | 2 | 6 | 6 | 0 | 0 | 1.3s |
| 2 | test_02_standard_fails | 3 | 12 | 8 | 2 | 1 | 33s |
| 3 | test_03_recurrence | 4 | 12 | 12 | 0 | 0 | 1.5s |
| 4 | test_04_spectral_gains | 5 | 8 | 8 | 0 | 0 | 28s |
| 5 | test_05_inexactness | 6 | 3 | 3 | 0 | 0 | 1.7s |
| 6 | test_06_voronoi | 7 | 7 | 5 | 2 | 0 | 3.1s |
| 7 | test_07_curvature | §4/§6 | 5 | 5 | 0 | 0 | 0.5s |
| **Total** | | | **53** | **47** | **4** | **1** | **~69s** |
