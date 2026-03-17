# Paper A1 — Argument Skeleton (v2)

**Title:** Exactness-preserving discrete de Rham complexes under Bloch-periodic boundary conditions
**One sentence:** Standard DEC breaks exact sequences under Bloch BCs; we fix it via recurrence and characterize the spectral consequences.

---

## Argument flow

### §1. Introduction (2 pages)

**Key sentence (must appear in intro):**
> "No construction is known to us that produces an exact discrete de Rham
> complex under Bloch-periodic boundary conditions on unstructured
> polyhedral meshes."

(Hedged: "no construction is known to us" not "no exact complex exists."
Backed by: Hirani 2003, Desbrun 2005, Schulz 2018, Mönkölä 2023 — none
address exactness under Bloch BCs.)

**Opening:** DEC is widely used for Maxwell on periodic structures.
Several references claim DEC avoids spurious modes (Schulz 2018, Mönkölä 2023).
FEEC (Arnold-Falk-Winther 2006, 2010) established that exact subcomplexes are
necessary for stable Hodge Laplacian discretization.

**The gap:** Under Bloch-periodic BCs on unstructured polyhedral meshes,
d₁(k)d₀(k) ≠ 0 for generic k. The exact sequence breaks. This is a
structural failure, not an implementation bug.

**Our contribution:**
1. We prove that standard DEC breaks exactness under Bloch BCs (Proposition 1)
2. We construct an exact d₁(k) via face-boundary recurrence (Theorem 1)
3. We prove uniqueness, canonical K, and kernel = V (Künneth)
4. We show that inexactness causes rank pollution, Hodge hybridization,
   and eigenvalue redistribution (§6)
5. As application, we identify Voronoi tessellations as geometries where
   the exact complex yields optimal dispersion (§7)

---

### §2. Setup (2 pages)

**Defines:**
- Periodic cell complex (V, E, F, C) on T³
- Incidence matrices d₀, d₁ (topological)
- Hodge stars ⋆₀, ⋆₁, ⋆₂ (metric, diagonal)
- Bloch-twisted operators: d₀(k)[e,v] = d₀[e,v]·exp(ik·n_e·L)
- Curl-curl operator: K(k) = d₁(k)†⋆₂ d₁(k)
- Test structures (table): Kelvin, C15, WP, SC, random Voronoi

No claims. Pure setup.

---

### §3. The exactness problem (3 pages)

**Proposition 1:** Standard Bloch curl d₁_std breaks exactness.
- Statement: If face f has ≥2 boundary edges with distinct lattice shifts
  sharing a vertex, then d₁_std(k)d₀(k) ≠ 0 for generic k.
- Proof: phase cancellation argument (3 lines)

**Corollary 1:** No per-edge phase assignment can fix this.
- Proof: intra-face contradiction (same-shift edges forced to different phases)

**Numerical evidence:**
- ‖d₁_std d₀‖ ~ O(1) on all test structures (table)
- n_spur = 3-16 spurious modes depending on structure/direction (table)

---

### §4. Exactness-preserving construction (4 pages)

**Lemma 1:** Flat holonomy. H_f = ∏ phases around face = 1.
- Proof: contractibility of face in R³ → zero net lattice translation

**Theorem 1 — MAIN RESULT:** Recurrence construction.
- Statement: For each face f with boundary edges e₀,...,e_{n-1}:
  set φ₀=1, then φᵢ = -σᵢ₋₁φᵢ₋₁ d₀[eᵢ₋₁,vᵢ]/(σᵢ d₀[eᵢ,vᵢ]).
  Then d₁(k)d₀(k) = 0 for all k.
- Proof: substitution into vertex equations + Lemma 1 for closure
- Cost: O(Σ n_f), same as standard d₁

**Corollary 2:** Uniqueness up to per-face gauge λ_f.

**Proposition 2:** Canonical K. Any two constructions give same K.

**Proposition 3:** Kernel dimension = |V| for generic k.
- Proof: 3 steps (injectivity of d₀, exactness, Künneth on T³).

**Remark:** Bloch cohomology β(Γ) = (1,3,3,1), β(k≠0) = (0,0,0,0).

**Formal result count for §3-§4:**
3 Propositions, 1 Theorem, 1 Lemma, 2 Corollaries. All with proofs.

---

### §5. Spectral consequences (3 pages)

Numerical results, not formal claims. Presented as "Results" not "Theorems."

**Result 1:** Spectral pollution eliminated.
- Exact: n_zero = V at all k, all structures, all directions.
- Standard: n_zero = V − n_spur. Table: 5 structures × 3 directions.

**Result 2:** Hodge splitting restored.
- Exact: gradient overlap < 10⁻¹² on all physical modes.
- Standard: 85-97% of modes mixed. Table: 3 structures.

**Result 3:** Dispersion convergence.
- Exact: c² → 1 as O(1/N²), p=2.00, R²=1.0000.
- Standard: c² oscillates, no convergence. Figure: log-log.

**Result 4:** Band structure.
- Figure: Γ-X-M-R-Γ exact vs standard.

---

### §6. Consequences of inexactness (2 pages, condensed)

Framing: "What goes wrong WITHOUT exactness — structural explanation."
NOT a separate theory. Consequences of §3 failure, explaining §5 numerics.

**3 results only** (strongest from W3):

**Result 5:** Rank pollution = rank excess (algebraic identity).
- rank(d₁_std) − rank(d₁_exact) = n_spur. Not empirical.

**Result 6:** Hodge hybridization converges to random baseline.
- Standard: gradient overlap → nV/nE ≈ 0.51. Universal (3 structures).
- Deep non-perturbative: not a small correction.

**Result 7:** Trace conservation.
- Σ(eigenvalue shifts) + Σ(spurious eigenvalues) = 0.
- Total trace is topological invariant.

**Dropped** (to avoid stealing focus):
- Level spacing Poisson→GUE → separate paper or appendix
- IPR surface/bulk scaling → separate paper
- Perturbation regimes → separate paper

---

### §7. Application: Voronoi tessellations (2-3 pages)

Framing: NOT "We prove Voronoi optimality."
YES: "This illustrates the effect of exactness in an optimal metric setting."

**Proposition 4:** G = Vol·I on any periodic Voronoi tessellation.
- Proof: divergence theorem on dual cells. (1 page, self-contained.)
- Explicit link: "Proposition 4 provides the metric ingredient; combined with
  Theorem 1 (exactness), this yields exact isotropic dispersion."

**Result 8:** Combined with Theorem 1 → ω² = |k|² + O(|k|⁴).
- Numerical verification on 3 cubic + 5 random Voronoi.
- Standard DEC: c² = 0.52-1.68, direction-dependent.
- Note: "Full spectral analysis of the Voronoi case, including
  the converse and moment conservation, will appear in [A2]."

---

### §8. Numerical experiments (2 pages)

Consolidated figures and tables from §5 + §7.
- Fig 1: Mesh wireframes (Kelvin + random Voronoi)
- Fig 2: Band structure comparison (exact vs standard)
- Fig 3: Convergence log-log
- Table: universality across 5+ structures

---

### §9. Discussion (1 page)

- FEEC connection: exact subcomplex on polyhedral mesh with diagonal Hodge stars.
  Not simplicial, not Whitney forms — complementary to existing FEEC constructions.
- What Voronoi adds: metric isotropy sufficient for c²=1. Necessity → [A2].
- Open: non-flat connections (magnetic flux), non-abelian gauge, non-T³ topology.
- 2D validation with photonic crystal benchmarks → [Paper B].

---

## Final formal result count

| Type | Count | Where |
|------|-------|-------|
| Theorem | 1 | §4 (recurrence — THE result) |
| Proposition | 4 | §3 (failure, no fix), §4 (canonical K, kernel), §7 (G=Vol·I) |
| Lemma | 1 | §4 (holonomy) |
| Corollary | 2 | §3 (no per-edge fix), §4 (uniqueness) |
| Numerical Results | 8 | §5 (4), §6 (3), §7 (1) |
| **Total formal** | **8** | |

8 formal results, not 17. Much cleaner for reviewer.

---

## Claim dependency graph

```
                    ┌─── Cor 2 (unique)
                    ├─── Prop 2 (canonical K)
Lem 1 ──→ THM 1 ───┤
                    ├─── Prop 3 (kernel = V)
                    ├─── Results 1-4 (spectral gains)
                    └─── Result 8 (c²=1, with Prop 4)

Prop 1 ──→ Cor 1 (no per-edge fix)
       ──→ Results 5-7 (inexactness consequences)

Prop 4 (G=Vol·I, independent) ──→ Result 8 (c²=1)
```

Two roots: Theorem 1 (positive) and Proposition 1 (negative).
One independent geometric result: Proposition 4 (Voronoi).
Everything else follows.

---

## What reviewer sees

**Quick scan (2 min):**
- Title: exact de Rham under Bloch — clear, specific
- Abstract: problem → solution → consequences
- 1 theorem, 4 propositions — not overloaded
- Voronoi in §7 as application, not core

**Detailed read (1 hour):**
- §3: "I didn't know standard DEC fails here" → interest
- §4: proofs are algebraic, short, verifiable → confidence
- §5: clean numerics, exact vs standard contrast → convinced
- §6: "ah, so THAT's why standard is so bad" → understanding
- §7: "nice that Voronoi works, natural example" → completeness
- §9: connects to FEEC → "this person knows the field"

**Verdict:** "Clear contribution. One main theorem. Well-supported."

---

## New findings from test suite (17 Mar 2026)

These emerged during testing and may inform paper text:

1. **n_spur ordering [100]<[110]<[111] universal** on all 3 cubic structures.
   Kelvin (6,12,14), C15 (9,15,16), WP (3,5,7). Combinatorial origin: more
   boundary-crossing edges along higher-symmetry directions. Remark for §3.

2. **n_spur extensive on random Voronoi**: deficit 15 (V=331), 19 (V=347).
   Consistent with W3 R13 (pollution ∝ surface area). §3 material.

3. **Exactness purely topological**: ‖d₁d₀‖ = 7.86e-16 identically at ε=0%..50%
   mesh perturbation. Zero geometric dependence. Strong statement for §4.

4. **n_zero > V at BZ boundary (TRIM)**: n_zero=98 at k=k_BZ, rank(d₀)=95.
   Expected from Künneth (H¹≠0 when e^{ikL}=1). Documented, not a bug. §4 remark.

5. **Standard c² oscillates with N**: 1.55→2.48→1.93→1.03 (Kelvin N=2..5).
   Non-convergent. c²=1.03 at N=5 is accidental. §5 material.

6. **c²_std on standard DEC is SPURIOUS, not acoustic**: first non-zero eigenvalue
   on standard DEC is an expelled gauge mode (88-94% gradient). Standard cannot
   even identify the acoustic mode. More dramatic than just c²≠1. §5 material.

7. **Trace conservation at all k**: tr(M⁻¹K) rel_diff = 1.78e-16 across 5
   k-directions. Topological invariant, not coincidence. §6 material.

8. **Both ⋆₁ and ⋆₂ matter for c²=1**: perturbing ⋆₂ gives Δc²=0.024,
   perturbing ⋆₁ gives Δc²=0.012. Neither alone sufficient. §7 material.

9. **‖d₁_std d₀‖ ~ O(k⁰·⁹⁶)** at small k (3 points), O(k⁰·⁹⁰) on full range.
   Confirms O(k) scaling with higher-order correction. §3 material.
