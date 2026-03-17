# Tests Architecture — Paper A1

Tests organized by argument flow. Each test group supports a paper section.
Source indicates where existing code lives (to be ported).

---

## GROUP 1: Setup validation (§2)

Verify test structures are correctly built and DEC operators are well-formed.

| Test | What | Source | Priority |
|------|------|--------|----------|
| T1.1 | d₀ has correct shape (E×V), entries ±1 | bloch 1_test_core | Must |
| T1.2 | d₁_topological: d₁d₀ = 0 at k=0 | bloch 1_test_core | Must |
| T1.3 | ⋆₁, ⋆₂ positive diagonal | bloch 5_test_r16 (Part 0) | Must |
| T1.4 | Structures: V,E,F,C counts match expected | bloch 1_test_core | Must |
| T1.5 | Scalar Laplacian K₀ = d₀†⋆₁d₀ has no spurious zeros | bloch 4_test_structure | Nice |

**Structures:** Kelvin N=2, C15 N=1, WP N=1, SC N=3, random Voronoi (3 seeds).

---

## GROUP 2: The problem — standard fails (§3)

Demonstrate that standard Bloch d₁ breaks exactness.

| Test | What | Claim | Source | Priority |
|------|------|-------|--------|----------|
| T2.1 | ‖d₁_std(k) d₀(k)‖ > 1 on all structures | Prop 1 | bloch 1_test_core | Must |
| T2.2 | ‖d₁_std d₀‖ direction-dependent (axis < face < body) | Prop 1 | bloch 3_test_robustness | Must |
| T2.3 | n_spur = V − n_zero_std > 0 on all structures | Prop 1 | bloch 1_test_core | Must |
| T2.4 | Intra-face contradictions count (Kelvin: 37/112, SC: 45/81) | Cor 1 | bloch 8_test_cor1 | Must |
| T2.5 | Per-face constraint has 1D null space (uniqueness premise) | Cor 1 | bloch 8_test_cor1 (Part 0) | Must |
| T2.6 | ‖d₁_std d₀‖ ~ O(k) scaling | — | bloch 5_test_r16 (Part 4) | Nice |
| T2.7 | **NEG:** Standard ker(K) ≠ V (ker deficit = n_spur) | Prop 3 contrast | bloch 7_test_ker | Must |
| T2.8 | **NEG:** Standard cohomology undefined (d₁d₀≠0 → not a complex) | §4.6 contrast | bloch 13_test_cohomology | Must |
| T2.9 | **CONTROL:** At k=0, exact = standard (no Bloch → no problem) | Prop 1 scope | bloch 5_test_r16 (Part 0) | Must |

---

## GROUP 3: The fix — recurrence construction (§4)

Verify the construction and its algebraic properties.

| Test | What | Claim | Source | Priority |
|------|------|-------|--------|----------|
| T3.1 | ‖d₁_exact(k) d₀(k)‖ < 10⁻¹⁴ on all structures, all k | Thm 1 | bloch 1_test_core | Must |
| T3.2 | Holonomy H_f = 1 on all faces, all structures | Lem 1 | **NEW** (extract from gauge_bloch.py) | Must |
| T3.3 | K canonical under face permutations (6 types) | Prop 2 | bloch 3_test_robustness | Must |
| T3.4 | K invariant under gauge transform e^{iθ(v)} | Prop 2 | bloch 3_test_robustness | Must |
| T3.5 | n_zero = V at all k, all structures | Prop 3 | bloch 1_test_core | Must |
| T3.6 | rank chain: rank(d₀) = V, rank(d₁) = E−V, ker(d₁) = im(d₀) | Prop 3 | bloch 7_test_ker | Must |
| T3.7 | Exactness preserved under mesh perturbation ε=0..20% | Thm 1 | bloch 20_test_perturbation | Must |
| T3.8 | d₂d₁ = 0 (full complex extends to level 2→3) | §4.6 | bloch 9_test_d2 | Must |
| T3.9 | β(Γ) = (1,3,3,1), β(k≠0) = (0,0,0,0) | Claim 8 | bloch 13_test_cohomology | Must |
| T3.10 | BZ boundary behavior (rank drop at k=k_BZ, exactness maintained) | — | bloch 3_test_robustness | Nice |

---

## GROUP 4: What you gain — spectral consequences (§5)

Show the concrete benefits of exactness.

| Test | What | Claim | Source | Priority |
|------|------|-------|--------|----------|
| T4.1 | Hodge splitting: exact 0 mixed modes, standard 169-187/192 mixed | Claim 10 | bloch 1_test_core | Must |
| T4.2 | Gradient overlap < 10⁻¹² on physical modes (exact) | Claim 10 | bloch 7_test_ker | Must |
| T4.3 | Gradient overlap ~ 0.51 on physical modes (standard) | Claim 10 | bloch 7_test_ker | Must |
| T4.4 | c² convergence: exact O(1/N²) p=2.00, standard non-convergent | Claim 11 | bloch 2_test_convergence | Must |
| T4.5 | c² convergence on BOTH Kelvin and SC (two families) | Claim 11 | bloch 2_test_convergence | Must |
| T4.6 | Band structure Γ-X-M-R-Γ: exact clean, standard contaminated | Claim 12 | bloch 6_make_figures | Must |
| T4.7 | Universality: same results on 5+ structures (table) | — | bloch 1_test_core + 4_test_structure | Must |
| T4.8 | Random Voronoi: 10/10 seeds pass exactness | — | bloch 4_test_structure | Must |

---

## GROUP 5: Consequences of inexactness (§6, condensed)

3 results only — explains §5 numerics, doesn't steal focus from Theorem 1.

| Test | What | Result | Source | Priority |
|------|------|--------|--------|----------|
| T5.1 | rank(d₁_std) − rank(d₁_exact) = n_spur (algebraic identity) | R5 | W3 file 3 (R7-R8) | Must |
| T5.2 | Hybridization: gradient overlap → nV/nE ≈ 0.51 (3 structures) | R6 | W3 file 5 (R12) | Must |
| T5.3 | Trace conservation: Σ(shift) + Σ(spur) = 0 | R7 | W3 file 6 (R15b) | Must |

Dropped from paper (saved for future):
- Level spacing Poisson→GUE (W3 R15a)
- IPR surface/bulk scaling (W3 R13)
- Non-perturbative regime (W3 R14)

---

## GROUP 6: Voronoi optimality (§7)

Condensed from voronoi_maxwell. Only essential tests.

| Test | What | Claim | Source | Priority |
|------|------|-------|--------|----------|
| T6.1 | G = Vol·I on Kelvin, C15, WP (ε_G < 10⁻¹⁵) | Claim 16 | voronoi 1_test_metric | Must |
| T6.2 | G = Vol·I on 5 random Voronoi seeds | Claim 16 | voronoi 1_test_metric | Must |
| T6.3 | H = Vol·I on same structures | Claim 16 | voronoi 1_test_metric | Must |
| T6.4 | c² = 1 on Kelvin, C15, WP (exact DEC) | Claim 17 | voronoi 2_test_c_squared | Must |
| T6.5 | c² = 1 on 5 random Voronoi (exact DEC) | Claim 17 | voronoi 2_test_c_squared | Must |
| T6.6 | c²_std ≠ 1 (1.25-1.68 on cubic structures) | Claim 17 | voronoi 3_test_removing | Must |
| T6.7 | Removing Voronoi geometry → c ≠ 1 (perturbed stars) | — | voronoi 4_test_removing_geo | Nice |

---

## Summary

| Group | §  | Tests | Must | Nice | New code needed? |
|-------|-----|-------|------|------|-----------------|
| 1. Setup | 2 | 5 | 4 | 1 | No (port) |
| 2. Problem | 3 | 6 | 5 | 1 | No (port) |
| 3. Fix | 4 | 10 | 9 | 1 | No (port) |
| 4. Gains | 5 | 8 | 8 | 0 | No (port) |
| 5. Inexactness | 6 | 3 | 3 | 0 | **Yes** (port W3) |
| 6. Voronoi | 7 | 7 | 6 | 1 | No (port) |
| **Total** | | **39** | **35** | **4** | |

---

## Source mapping

| Source project | Files to port | Tests used |
|---------------|--------------|------------|
| st_bloch_exactness/tests/ | 1,2,3,4,5,6,7,8,9,13,20 | Groups 1-4 |
| st_bloch_exactness/internal/w_3/ | 3,5,6 | Group 5 (new) |
| st_voronoi_maxwell/tests/ | 1,2,3,4 | Group 6 |
| st_bloch_exactness/src/ | gauge_bloch.py, bloch.py | All (shared code) |
| st_base/src/ | core_math/*, physics/* | All (shared code) |

---

## File organization (proposed)

```
tests/
├── 01_test_setup.py           ← Group 1 (T1.1-T1.5)
├── 02_test_standard_fails.py  ← Group 2 (T2.1-T2.6)
├── 03_test_recurrence.py      ← Group 3 (T3.1-T3.10)
├── 04_test_spectral_gains.py  ← Group 4 (T4.1-T4.8)
├── 05_test_inexactness.py     ← Group 5 (T5.1-T5.3) **NEW**
├── 06_test_voronoi.py         ← Group 6 (T6.1-T6.7)
└── 07_make_figures.py         ← Figures for paper
```

---

## What's NOT tested here (parked for Paper A2 or B)

| Content | Why not |
|---------|---------|
| Dielectric ε (tests 10,12,15,16,17) | Physics application, not core math |
| MPB comparison | Needs external dependency |
| h-refinement (test 11) | Covered by k-convergence (test 2) |
| Spectral pairing (test 18) | Nice but not core claim |
| Hodge decomposition SVD (test 19) | Nice but overlap < 10⁻¹² already shows this |
| Acoustic deformation (test 14) | Physics interpretation |
| 2D results (W24) | Paper B |
| Moment conservation (voronoi test 5) | Paper A2 |
| Converse c²=1 → G=Vol·I | Paper A2 |
| Error factorization (voronoi test 4) | Paper A2 |
