# Gap Analysis — Tests Not Ported

Systematic comparison: st_bloch_exactness tests (20 files) + st_voronoi_maxwell
tests (7 files) vs st_exact_deRham tests (6 files).

**Question:** Are there tests in source projects that cover Paper A1 claims
but are missing from our test suite?

---

## Source tests → our coverage

### PORTED (covered by our Groups 1-6)

| Source | Our test | Status |
|--------|----------|--------|
| bloch 1_test_core: exactness, pollution, c², leakage | T2.1, T2.3, T4.2, T4.3, T4.7 | ✓ |
| bloch 2_test_convergence: O(1/N²) | T4.4, T4.5 | ✓ |
| bloch 3_test_robustness: directions, perturbation, gauge, BZ boundary | T2.2, T3.3, T3.4, T3.7, T3.10 | ✓ |
| bloch 4_test_structure: random Voronoi, n_scaling, scalar K₀ | T4.8, T1.5 | ✓ |
| bloch 7_test_ker_theorem: rank chain, standard failure | T3.5, T3.6, T2.7 | ✓ |
| bloch 8_test_cor1_proof: uniqueness, contradictions | T2.4, T2.5 | ✓ |
| bloch 9_test_d2_exactness: d₂d₁ | T3.8 (PARTIAL) | ✓ partial |
| bloch 13_test_bloch_cohomology: Betti numbers | T3.9 (PARTIAL) | ✓ partial |
| bloch 20_test_perturbation_stability: ε=0..20% | T3.7 | ✓ |
| W3 file 3: rank/spectral flow | T5.1 | ✓ |
| W3 file 5: hybridization | T5.2 | ✓ |
| W3 file 6: trace/moments | T5.3 | ✓ |
| voronoi 1_test_metric: G=Vol·I, H=Vol·I | T6.1, T6.2, T6.3 | ✓ |
| voronoi 2_test_c_squared: c²=1 | T6.4, T6.5 | ✓ |
| voronoi 3_test_removing_exactness: c²_std≠1 | T6.6 | ✓ |
| voronoi 4_test_removing_geometry: perturbed stars | T6.7 | ✓ |

### NOT PORTED — intentionally excluded (Paper A2 or B material)

| Source | What | Why excluded |
|--------|------|-------------|
| bloch 10_test_dielectric | ε-dependent tests | Physics application, not core math |
| bloch 11_test_h_refinement | h-refinement on Voronoi | Convergence covered by T4.4/T4.5 |
| bloch 12_test_mpb_comparison | MPB benchmark | External dependency, not core |
| bloch 14_test_acoustic_deformation | Helmholtz splitting | Physics interpretation |
| bloch 15_test_hodge_star_interface | Interface averaging | Paper A2 material |
| bloch 16_test_voronoi_dielectric | Dielectric on random Voronoi | Physics application |
| bloch 17_test_face_normal_statistics | cos²(n̂,k̂) statistics | Paper A2 material |
| bloch 18_test_spectral_pairing | ±λ pairing (Dirac) | Nice but not core claim |
| bloch 19_test_hodge_decomposition | SVD Hodge decomposition | Covered by T4.1 gradient overlap |
| voronoi 5_test_spectral: moments | Moment conservation | Paper A2 material |
| voronoi 6_test_appendix_a: naive PT | Perturbation theory | Paper A2 material |
| voronoi 7_test_appendix_b: std anatomy | Mode mixing | Covered by T4.3 |

### NOT PORTED — potential gaps

| Source | What | Risk | Action |
|--------|------|------|--------|
| bloch 3: `test_random_directions()` 20 random dirs | We test 3 dirs (T2.2). 20 is more robust. | Low — 3 dirs already show direction-dependence. | **Note, don't port** |
| bloch 3: `test_operator_norm()` ‖K_std−K_exact‖_F vs k | Shows operator distance grows with k. | Low — not a paper claim. | **Note, don't port** |
| bloch 3: `test_asymptotics()` c² from k=10⁻⁴ to 0.20 | Full k-scan of c² stability. | Low — T4.4 covers convergence. | **Note, don't port** |
| bloch 4: `test_voronoi_scaling()` n=50..400 | c² vs Voronoi size. | Medium — shows c²→1 improves with mesh quality. | **Consider for T6** |
| bloch 4: `test_n_scaling()` n_spur grows with N | Pollution is extensive. | Medium — paper mentions this in §3. | **Consider** |
| bloch 5: `test_r16_oscillations` Parts 3-6 | k-scaling of spurious eigs, Hodge defect | Low — covered by T2.6 (slope). | **Note** |
| bloch 9: full d₂ Bloch recurrence | d₂(k)d₁(k)=0 exact | **HIGH** — T3.8 is PARTIAL. | **TODO (tracked)** |
| bloch 13: `test_bz_scan()` 36 k-points | Fine BZ scan of Betti | Low — T3.10 covers BZ boundary. | **Note** |
| bloch 13: `test_rank_structure()` | Rank identities at Γ vs k≠0 | Low — covered by T3.6. | **Note** |
| bloch 20: `test_2_betti_stable()` | Betti under perturbation | Medium — we test exactness under pert (T3.7) but not Betti. | **Consider** |
| bloch 20: `test_3_hodge_stable()` | Hodge decomp under perturbation | Medium — same as above. | **Consider** |

---

## Verdict

**No critical gaps.** All Paper A1 claims (Prop 1, Cor 1, Thm 1, Prop 2, Prop 3,
Results 1-8) are covered by existing tests.

**Two medium-priority gaps:**

1. **n_spur extensive** (bloch 4 `test_n_scaling`): n_spur grows with N (6→17→29).
   Paper §3 mentions pollution is extensive. We show it's direction-dependent (T2.2)
   and structure-dependent (T2.1) but not N-dependent.
   → **Add to T2 or note as §3 material from existing data**

2. **Betti stable under perturbation** (bloch 20): We test exactness under pert (T3.7)
   but not that Betti numbers are stable. Since exactness implies correct Betti
   (by construction), this is logically redundant but numerically valuable.
   → **Low priority — exactness stability implies Betti stability**

**One high-priority gap (already tracked):**

3. **d₂ Bloch recurrence** (T3.8 PARTIAL): Full complex d₂(k)d₁(k)=0.
   Already in tracker as TODO ~2-3h.
