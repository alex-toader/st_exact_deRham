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
