# st_exact_deRham — Paper A1 Tracker

**Date:** 17 Mar 2026
**Paper:** "Exactness-preserving discrete de Rham complexes under Bloch-periodic boundary conditions"
**Venue:** SIAM Journal on Numerical Analysis
**Target submit:** after 13 Apr 2026 (JCP transfer expiry)

---

## Structure

| § | Title | Pages | Source | Status |
|---|-------|-------|--------|--------|
| 1 | Introduction | 2 | Rewrite | TODO |
| 2 | Periodic cell complexes and DEC operators | 2 | JCP §2 | Port |
| 3 | The exactness problem | 3 | JCP §3 | Port |
| 4 | Exactness-preserving construction | 4 | JCP §4 | Port |
| 5 | Spectral consequences | 3 | JCP §5 | Condense |
| 6 | Consequences of inexactness | 3 | W3 (new) | TODO |
| 7 | Voronoi as optimal geometry | 3 | voronoi_maxwell | Condense |
| 8 | Numerical experiments | 2 | JCP figs/tables | Regenerate |
| 9 | Discussion | 1 | — | TODO |

## Formal results

| ID | Statement | Proof | Source |
|----|-----------|-------|--------|
| Prop 1 | Standard d₁ breaks exactness under Bloch | ✓ | JCP |
| Cor 1 | No per-edge phase fix exists | ✓ | JCP |
| Lem 1 | Flat holonomy H_f = 1 | ✓ | JCP |
| **Thm 1** | **Recurrence → d₁d₀=0** | ✓ | JCP |
| Cor 2 | Uniqueness up to per-face gauge | ✓ | JCP |
| Prop 2 | Canonical K (gauge-invariant) | ✓ | JCP |
| Prop 3 | ker K = V (Künneth) | ✓ | JCP |
| Thm (§7) | G=Vol·I on periodic Voronoi | ✓ | voronoi_maxwell |

## Key numerics needed

| Test | What | Source |
|------|------|--------|
| T1 | ‖d₁d₀‖ on 5 structures (exact vs standard) | JCP table 3 |
| T2 | O(h²) convergence log-log | JCP table 4, fig 3 |
| T3 | Hodge splitting overlap (exact < 10⁻¹⁵, standard 0.51) | JCP table 5 |
| T4 | Band structure (exact vs standard) | JCP fig 2 |
| T5 | Rank pollution = rank excess | W3 (new) |
| T6 | Hybridization at random baseline | W3 (new) |
| T7 | Trace conservation | W3 (new) |

## Bibliography plan

Must cite (FEEC): Arnold-Falk-Winther 2006, 2010; Christiansen 2008
Must cite (DEC): Hirani 2003; Desbrun-Leok-Marsden 2005
Must cite (polyhedral DEC): Gillette 2011; Mönkölä 2023; Schulz 2018
Must cite (spectral pollution): Boffi 2010; Caorsi 2001
Must cite (Maxwell DEC): Bossavit 1998
Target: 25-35 refs

## Open issues

- [ ] Nédélec comparison table (L²/H(curl) errors on standard problem). SIAM will likely request.
      Not blocking for first draft but should be in submitted version.
      Needs: FEniCS setup, Nédélec P1 on same T³ domain, eigenvalue comparison.
- [x] Abstract positioning: "standard DEC fails under Bloch" in first sentence. Draft 2 done.
- [ ] Prop 1 proof must be impeccable — reviewer will test "is this trivial or new?"
      Current proof is 3 lines. May need explicit example (face with 2 shifts) to convince.
- [ ] "No construction is known to us" — hedge version. Back with 4 citations that don't do it.
- [ ] §7 Voronoi link: explicit sentence connecting Prop 4 + Thm 1 → c²=1.
- [ ] Review note: HAS_VORONOI_DUAL mutable global — minor, idempotent, not blocking.
- [ ] Review note: T1.5 single k magnitude — T3.10 covers BZ boundary, not duplicating.
- [ ] Review note: T1.7 edge-in-2-faces — builders have internal asserts, not duplicating.
- [ ] Review note: K₀ isotropy identical on 3 dirs for cubic (Oh symmetry), varies on random.
      Scalar isotropy ≠ vector isotropy. G=Vol·I implies vector (1-form). Remark for §7.
- [x] DONE: n_spur ordering [100]<[110]<[111] universal — verified in T2.2 + T2.10.
      n_spur ~ N^2.29 ([100]), N^2.37 ([111]). Surface scaling. Remark for §3.
- [x] DONE: n_spur extensive on random Voronoi — verified in T2.1 (deficit 15, 19).
- [x] DONE: T2.11 per-edge optimization — min ‖d₁d₀‖ = 1.2564, ratio 10¹⁵ vs exact.
      Unique global minimum (10 starts converge to same). Cor 1 confirmed numerically.
- [x] **DONE:** d₂(k) Bloch recurrence implemented (src/physics/bloch_complex.py).
      Our version is improved vs st_base: adds minimum image convention (bug fix),
      holonomy check on non-BFS edges, edge safety in face_edge_map.
      st_base bloch_complex.py has minimum image bug in build_cell_face_incidence.
      T3.8 PASS: d₂d₁=0 at 10⁻¹⁵ on all 4 structures.
      T3.9 PASS: β(Γ)=(1,3,3,1) and β(k≠0)=(0,0,0,0) on all 4 structures.
      No more PARTIAL tests. Full cochain complex verified.
- [x] T3.7 confirmed topological: ‖d₁d₀‖ identical at ε=0%..50%. Pure topology.
- [x] T3.10: n_zero=98 at k_BZ explained — TRIM cohomology (H¹≠0 when e^{ikL}=1).
- [x] T4 threshold floor unified to 1e-10 (was 1e-14 in 3 places).
- [x] T4.4: c²_std=1.03 at N=5 documented as accidental (spread 1.45 confirms non-convergence).
- [ ] Review note: T4.9 (c²=1 random Voronoi) already covered by T6.5. No duplication.
- [ ] Review note: T4.10 (standard fails on random Voronoi) already covered by T2.1. No duplication.
- [ ] Tech debt: grad_overlap uses np.linalg.inv — replace with solve/Cholesky for large structures.
- [ ] Paper §6 material: T5.3 multi-k verification (1.78e-16 across 5 directions) confirms
      trace conservation is topological invariant. Use in paper: "holds at all k tested."

## LaTeX formatting notes

- [ ] Master theorem (i)-(iv): restructure as Theorem 1 (existence) + Corollaries 1-3
      for uniqueness, canonical K, kernel. Standard SIAM format. Do in .tex, not .md.
- [ ] Proposition numbering: verify sequential in final .tex
- [ ] SM: format as separate .tex with shared preamble

## Progress

- [x] skeleton.md v2 (argument flow, 8 formal results)
- [x] tests_architecture.md (39 tests, 6 groups)
- [x] abstract draft 1
- [x] Master theorem added at start of §4 (unifies Lem 1 + Thm 1 + Cor 2 + Prop 2 + Prop 3)
- [ ] Strategic: reduce experimental by ~30% before submission (keep SC + Kelvin + 1 random,
      drop C15/WP redundancy in some sections). Reviewer: "I believe you after 2 meshes, not 5."
- [x] Strategic: Voronoi de-escalated to Structural property + SM. Feeds Paper A2.
- [ ] Paper A2 material ready: voronoi_maxwell (42 tests, 5 thms) + W22 (40 results) + lorentz (70 results).
      Venue: Discrete Comput. Geom. Content: G=Vol·I full proof, converse, moments, hierarchy.
- [ ] Strategic: promote Prop 5 (rank-pollution) earlier or make more central.
- [x] DONE: T3.11 iff completeness — nullity = nF on all 4 structures. Recurrence is ONLY solution.
- [x] DONE: T2.12 sensitivity — ‖d₁d₀‖/ε = 33.71 constant (CV=0.001). Exact is isolated point.
## Reviewer ideas (triaged 17 Mar 2026)

- [ ] A2: O(k⁴) dispersion coefficient a₄ — fit across structures, check if geometric. 30 min.
- [ ] A3: n_spur along BZ path Γ→X→M→R→Γ — continuous plot, would be strong Figure 3(c). 1h.
- [ ] B1: Obstruction class ∈ H²(T³, U(1)) — cohomological interpretation of ‖d₁_std d₀‖. Paper A2.
- [ ] B2: Generalization to T^d — recurrence works on any T^d. Theorem for Paper B intro.
- [ ] C1: Inverse design — optimize mesh for max gap within admissible perturbations. Paper A2/standalone.
- [ ] C2: Discrete Dirac operator — ±λ pairing from exact complex. Standalone paper.
- [ ] C3: Topological band gaps on 3D — killer app if exact DEC detects gaps standard misses. Paper B or later.

## Reviewer questions answered

- [x] T7 holonomy order: ordered vs unordered product identical to 10⁻¹⁶.
      U(1) is abelian — complex multiplication commutes. Reviewer concern invalid.
- [x] T7 np.angle wrapping: not a problem at frac ≤ 0.10 (max individual phase ~36°).
      Documented in code. Would need np.unwrap at large k.
- [x] T7 n_curv > n_spur analogy: softened from "field vs charge" to "partial cancellation
      observation." Not a theorem, just a numerical fact. Documented.
- [x] Künneth at TRIM: β=(0,0,0,0) at X,M,R (λ=-1≠1, acyclic). BZ_x=(2π/L,0,0)
      is equivalent to Γ (λ=1, β=(1,3,3,1), n_zero=98). Paper corrected: -1≠1.
- [x] Voronoi near-degeneracy: G=Vol·I holds (err=5.8e-12) even with 4 clustered
      points (min star1=0.000025). General position needed for BUILDER, not for identity.
- [x] T7 support inclusion: [100] ⊂ [110] exactly, [100]∩[111]=23/24, [110]∩[111]=39/40.
      No clean hierarchy. 52/112 faces have curvature on at least one direction.
      Interesting geometry but not a theorem. Noted, not in paper.

## Progress

- [ ] §1 Introduction drafted
- [ ] §2 ported from JCP
- [ ] §3 ported from JCP
- [ ] §4 ported from JCP
- [ ] §5 condensed from JCP
- [ ] §6 written (new)
- [ ] §7 condensed from voronoi_maxwell
- [ ] §8 figures regenerated
- [ ] §9 Discussion drafted
- [ ] Bibliography compiled
- [ ] Internal review
- [ ] Submit
