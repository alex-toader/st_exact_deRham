# Exactness-preserving discrete de Rham complexes under Bloch-periodic boundary conditions

Alexandru Toader (Independent researcher)

---

## Abstract

Discrete exterior calculus (DEC) is often described as free of spurious
modes in the curl–curl eigenvalue problem. We show that this property fails
under Bloch-periodic boundary conditions: the standard discrete curl satisfies
d₁(k)d₀(k) ≠ 0 for generic wave vector k on any unstructured polyhedral mesh
with faces whose boundary edges carry distinct lattice shifts. No per-edge
phase assignment can restore exactness. We construct an alternative d₁(k) via
a face-boundary recurrence that enforces d₁(k)d₀(k) = 0 for all k. The
construction is explicit, unique up to per-face gauge, and produces a canonical
curl–curl operator with gauge kernel of dimension |V|. Without exactness, the
pollution count equals the rank excess of d₁ and the Hodge decomposition
degrades to random-baseline hybridization. Experiments on five polyhedral
complexes confirm that exactness eliminates spectral pollution, restores the
Hodge splitting, and yields second-order dispersion convergence. On periodic
Voronoi tessellations, exactness combined with metric isotropy yields
ω² = |k|² + O(|k|⁴), verified numerically and supported by analytical arguments.

---

## 1. Introduction

The discrete exterior calculus (DEC) provides a coordinate-free discretization
of differential forms on cell complexes [Hirani 2003, Desbrun 2005]. For
Maxwell's equations on periodic structures, DEC represents the electric field
as a 1-cochain (edge values) and the magnetic flux as a 2-cochain (face values),
with discrete gradient d₀ and curl d₁ acting between cochain spaces. The
curl–curl eigenvalue problem

    K a = ω² M a,    K = d₁† ⋆₂ d₁,    M = ⋆₁

yields electromagnetic eigenfrequencies when the discrete de Rham sequence

    C⁰ →^{d₀} C¹ →^{d₁} C²

is exact, i.e., d₁ d₀ = 0. Boffi [2010] established that exactness of the
discrete complex is necessary and sufficient for spurious-free Maxwell
eigenvalue computations in the finite element context. The argument extends
to DEC: exactness d₁d₀ = 0 implies ker K = im d₀ (gauge completeness), and
the positive-definite Hodge star ⋆₂ ensures discrete compactness of the
physical subspace.

To our knowledge, no construction has been reported that produces an exact
discrete de Rham complex under Bloch-periodic boundary conditions on
unstructured polyhedral meshes. Unlike standard periodic identification,
Bloch conditions introduce a nontrivial rank-1 local system parametrized by
k, so exactness is no longer purely combinatorial — it depends on the
compatibility of phases across face boundaries.

On structured grids (Yee lattice), exactness under Bloch BCs holds automatically
by the tensor-product topology. DEC-based Maxwell computations use either
structured grids [Teixeira 1999, Schulz 2018, Mönkölä 2023] or non-periodic
domains [Chen 2017]. The discrete de Rham (DDR) method [Di Pietro 2020] achieves
exactness on polyhedral meshes but has not been formulated with Bloch periodicity.
The FEEC framework [Arnold-Falk-Winther 2006, 2010] provides the theoretical
foundation for structure-preserving discretizations but treats simplicial elements,
not polyhedral DEC.

Under Bloch-periodic boundary conditions on unstructured polyhedral meshes, the
standard construction of d₁(k) breaks exactness: d₁(k)d₀(k) ≠ 0 for generic k.
The failure has a combinatorial origin: faces whose boundary edges carry distinct
lattice shifts create intra-face phase contradictions that no per-edge assignment
can resolve.

**Contribution.** We present:
1. A proof that the standard Bloch extension breaks exactness on unstructured
   meshes, and that no per-edge phase fix exists (§3).
2. An explicit recurrence construction of an exactness-preserving d₁(k) that
   derives curl phases from those in d₀(k) (§4, Theorem 1).
3. Proofs of uniqueness, canonical K, and kernel dimension = |V| via Künneth
   formula on T³ with local coefficients (§4).
4. Three structural consequences of inexactness: rank pollution = rank excess
   (algebraic identity), Hodge hybridization at random baseline, and trace
   conservation (§6).
5. Identification of Voronoi tessellations as geometries where the exact
   complex yields optimal dispersion ω² = |k|² + O(|k|⁴) (§7).
6. Verification on five polyhedral complexes including random Voronoi meshes
   with no symmetry (§5, §7).

**Viewpoint.** The core of the construction is the observation that Bloch
boundary conditions define a rank-1 local system L_k on T³, and the discrete
operators must form a cochain complex with coefficients in L_k. The standard
DEC approach assigns Bloch phases independently per edge — this is a discrete
connection with nonzero curvature. Our recurrence derives all phases from
d₀(k), producing a flat connection that preserves exactness. The paper can
thus be read as constructing a flat discrete connection on the Bloch local
system such that the twisted cochain complex remains exact.

**Relation to FEEC.** In the FEEC framework [Arnold-Falk-Winther 2006, 2010],
exactness is preserved via commuting projections on simplicial subcomplexes.
Bloch-periodic twisting introduces a local system structure not explicitly
addressed in existing FEEC constructions — the standard FEEC machinery handles
periodicity by mesh identification but does not explicitly address the phase consistency
required by nontrivial Bloch parameters. Our construction fills this gap for
polyhedral DEC with diagonal Hodge stars. We do not address the full FEEC
program (commuting projections, approximation theory); the contribution is
the algebraic exactness of the twisted discrete complex.

---

## 2. Periodic cell complexes and DEC operators

### 2.1 Cell complex on the torus

Let Λ = L₁Z × L₂Z × L₃Z be a lattice in R³ with fundamental domain
Ω = [0,L₁) × [0,L₂) × [0,L₃). A periodic cell complex (V, E, F, C) on the
flat torus T³ = R³/Λ consists of vertices V, oriented edges E, oriented
faces F, and 3-cells C.

Edges connecting vertices across the periodic boundary carry a lattice shift
vector n_e ∈ Z³:

    n_e = round((x_j − x_i) / L)

where division is componentwise. Interior edges have n_e = 0.

### 2.2 Incidence matrices

The topological gradient d₀ ∈ R^{|E|×|V|} and curl d₁ ∈ R^{|F|×|E|} encode
boundary relations:

    d₀[e, v] = +1 if v = head(e), −1 if v = tail(e), 0 otherwise.

    d₁[f, e] = +1 if e ∈ ∂f with same orientation, −1 if opposite, 0 if e ∉ ∂f.

These satisfy d₁ d₀ = 0 (every vertex appears in the boundary of a face an
even number of times with canceling signs).

### 2.3 Hodge stars

The diagonal Hodge star operators ⋆₁ ∈ R^{|E|×|E|} and ⋆₂ ∈ R^{|F|×|F|}
encode metric information from the dual complex:

    ⋆₁[e, e] = |σ*_e| / |σ_e|,    ⋆₂[f, f] = |σ*_f| / |σ_f|

where σ_e, σ_f are primal edge length and face area, and σ*_e, σ*_f are
their dual counterparts. For the SC cubic lattice with constant a,
the stars are uniform: ⋆₁ = a·I, ⋆₂ = a⁻¹·I.

### 2.4 Bloch-twisted operators

For wave vector k ∈ R³, the Bloch-twisted gradient acts on 0-cochains u by:

    (d₀(k) u)_e = exp(ik · n_e · L) · u_{head(e)} − u_{tail(e)}

When n_e = 0 (interior edge), this reduces to u_{head} − u_{tail} = (d₀ u)_e.
At k = 0, d₀(k) = d₀ (topological gradient).

The standard Bloch curl applies independent per-edge phases:

    d₁_std(k)[f, e] = d₁[f, e] · exp(ik · n_e · L)

This is the construction used in existing DEC-based Maxwell solvers. We show
in §3 that d₁_std(k) d₀(k) ≠ 0 for generic k on unstructured meshes.

### 2.5 Test structures

We verify the construction on five periodic polyhedral complexes spanning
three qualitatively different categories:

| Structure | Type | |V| | |E| | |F| | |C| |
|-----------|------|-----|-----|-----|-----|
| SC (N=3) | Regular cubic lattice | 27 | 81 | 81 | 27 |
| Kelvin (N=2) | BCC Voronoi tessellation | 96 | 192 | 112 | 16 |
| C15 (N=1) | Laves phase Voronoi | 136 | 272 | 160 | 24 |
| WP (N=1) | A15 Voronoi | 46 | 92 | 54 | 8 |
| Random Voronoi | 50 uniform points | 206–354 | 412–708 | 236–404 | 50 |

The Voronoi tessellations have 4 edges per vertex and 3 faces per edge
(Plateau structure). Random Voronoi meshes are constructed from 50 uniformly
distributed points on T³; the range reflects 10 different seeds.

[FIGURE 1: Wireframe views of Kelvin N=2 and random Voronoi (n=50).]

---

## 3. The exactness problem

### 3.1 Standard Bloch curl breaks exactness

**Proposition 1.** If face f has ≥2 boundary edges with distinct lattice shifts
n_a ≠ n_b sharing a vertex, then d₁_std(k)d₀(k) ≠ 0 for all k outside the
hyperplane {k : exp(ik·(n_a − n_b)·L) = 1}. In particular, the failure occurs
on an open dense set of k — it is structural, not numerical.

*Proof.* At the shared vertex v, the two contributions to (d₁d₀)[f,v] carry
phases exp(ik·n_a·L) and exp(ik·n_b·L). Cancellation requires
exp(ik·(n_a − n_b)·L) = 1, which defines a hyperplane in k-space. This holds
at k = 0 but fails generically when n_a ≠ n_b. □

**Numerical evidence:**

| Structure | ‖d₁_std d₀‖ | n_spur [100] | n_spur [111] |
|-----------|-------------|-------------|-------------|
| Kelvin N=2 | 7.54 | 6 | 14 |
| C15 N=1 | 7.80 | 9 | 16 |
| WP N=1 | 5.07 | 3 | 7 |
| SC N=3 | 6.21 | 5 | 13 |
| Random Voronoi | 10.5–11.7 | 15–19 | — |

Direction-dependent pollution: [100] < [110] < [111], universal on all structures.

### 3.2 No per-edge fix exists

**Corollary 1.** Let φ: E × R³ → C* be any function assigning a multiplicative
phase to each edge, depending only on the edge — i.e., independent of the
face context or position along the face boundary.
Define d̃₁(k)[f,e] = d₁[f,e] · φ(e,k). Then d̃₁(k) d₀(k) ≠ 0 for generic k
on any mesh with faces whose boundary edges carry distinct lattice shifts.

More precisely: no 1-cochain multiplicator can make d₁ exact. The obstruction
is that exactness requires face-dependent phases (each face determines its own
phase assignment via the recurrence), not edge-dependent phases.

*Proof.* By the uniqueness result of §4.3 (Corollary 2 below), any local d₁
with d₁d₀ = 0 has d₁[f,:] = λ_f · d₁_exact[f,:] — the phases
depend on the face, not just on the edge. A multiplicative 1-cochain
φ(e,k) assigns the same factor to edge e regardless of which face it belongs
to. But the recurrence gives ψ_b/ψ_a = exp(ik·S·L) where S ∈ Z³ is the net
lattice shift along the partial boundary path from position a to b. When two
edges of the same face have the same lattice shift n_e but different positions
in the face boundary, their exact phases differ by exp(ik·S·L) ≠ 1.
A 1-cochain multiplicator cannot produce this face-dependent variation.

*Explicit example.* On Kelvin N=2, face 1 has three shift-0 edges at positions
0, 2, 4. The recurrence gives ψ₂/ψ₀ = exp(−ik_x L) ≠ 1 for k_x ≠ 0. Three
edges with the same shift require three distinct phases — impossible for a
per-edge function. Out of 112 faces, 37 exhibit such contradictions. □

**Numerical demonstration.** Optimizing 192 free edge phases via L-BFGS-B
(10 random starts, all converge to unique global minimum):

    min_θ ‖d₁(θ) d₀‖ = 1.2564

while the recurrence construction (§4) achieves 7.86 × 10⁻¹⁶. Ratio: 10¹⁵.
Per-edge approach is fundamentally impossible, not merely suboptimal.

**Remark (extensive pollution).** The pollution count grows as n_spur ~ N^{2.3}
on Kelvin (N = 2, 3, 4), between surface (N²) and volume (N³) scaling.
This is consistent with spurious modes being surface states concentrated on
boundary-crossing edges (see also §6.1).

---

## 4. Exactness-preserving construction

### 4.1 Flat holonomy

**Lemma 1.** For any face f, the holonomy product H_f = ∏ phases around the
boundary = 1.

*Proof.* Choose a lift of face f to the fundamental domain of T³ in R³. The
lifted face is a contractible polygon whose boundary is a closed path in R³.
Each boundary edge contributes a lattice shift n_i to the holonomy product
H_f = exp(ik·(Σ n_i)·L). The lift is chosen so that consecutive edges
concatenate in R³ without jumps. Since the resulting path is closed in R³,
the total lattice translation is Σ n_i = 0, giving H_f = 1. (Faces crossing the periodic
boundary are handled by choosing any consistent lift; the result is
independent of the lift since different lifts differ by a full lattice
vector, which contributes exp(ik·N·L) = 1 to the holonomy.) □

Verified: max |H_f − 1| = 2.3 × 10⁻¹⁶ on all faces, all structures.

### 4.2 Recurrence construction

**Theorem 1 (Existence).** There exists a local operator d₁(k): C¹ → C² with
support on face boundaries satisfying d₁(k)d₀(k) = 0 for all k. For each
face f with ordered boundary edges e₀, ..., e_{n-1}, define d₁(k) by:

    Set φ₀ = 1.
    For i = 1, ..., n-1:
        φᵢ = -σᵢ₋₁ φᵢ₋₁ d₀(k)[eᵢ₋₁, vᵢ] / (σᵢ d₀(k)[eᵢ, vᵢ])
    Set d₁(k)[f, eᵢ] = σᵢ φᵢ.

*Proof.* The condition d₁(k) d₀(k) = 0 at face f decomposes into n vertex
equations. At vertex vᵢ for i = 1, ..., n−1:

    d₁[f, eᵢ₋₁] · d₀(k)[eᵢ₋₁, vᵢ] + d₁[f, eᵢ] · d₀(k)[eᵢ, vᵢ] = 0

Substituting d₁[f, eᵢ] = σᵢ φᵢ yields the recurrence. This determines
φ₁, ..., φ_{n-1} from φ₀ = 1. The remaining equation at v₀ involves
d₁[f, e_{n-1}] and d₁[f, e₀], and is satisfied if and only if the product
of all recurrence ratios around the face equals 1 — which is exactly the
holonomy H_f = 1 guaranteed by Lemma 1. (Equivalently: the n vertex equations
sum to zero by the topological identity d₁d₀ = 0 at k = 0, so any n−1
of them imply the remaining one.)

Cost: O(Σ n_f), the same as assembling the standard d₁. □

Conceptually, the construction enforces phase consistency along each face
boundary, effectively defining a flat discrete connection induced by d₀(k).
The flatness (Lemma 1) guarantees that this connection is well-defined
regardless of the traversal order.

**Verified:** max ‖d₁_exact d₀‖ = 9.53 × 10⁻¹⁶ across 4 structures × 8
k-points. Exactness identical under 0%–50% vertex perturbation (topological
property, independent of metric).

**Remark (isolated solution).** The exact construction is an isolated point
in phase space, not a basin of attraction. Perturbing the phases of d₁_exact
by ε (random direction) produces ‖d₁d₀‖ = O(ε) with constant ~34, linear
across 5 orders of magnitude (ε = 10⁻⁸ to 10⁻¹). Any deviation from the
exact phases destroys exactness proportionally.

### 4.3 Uniqueness

**Corollary 2 (Uniqueness).** Any local operator satisfying d₁d₀ = 0 with
support on ∂f differs from the recurrence construction by a scalar
λ_f ∈ C* per face: d̃₁[f,:] = λ_f · d₁[f,:]. There are no other solutions.

*Proof.* The exactness condition at face f is a linear recurrence of order 1
with n−1 equations in n unknowns. Once d₁[f, e₀] is chosen (the seed), all
other entries are determined. The solution space is therefore one-dimensional
per face. □

The free phase λ_f per face is a discrete analogue of gauge freedom: it does
not affect the curl–curl operator (Corollary 3) but changes the representation
of d₁.

**Verified (iff).** The full constraint matrix for d₁d₀ = 0 with fixed support
on ∂f has nullity exactly nF on all 4 test structures (Kelvin: 112, C15: 160,
WP: 54, SC: 81). There are no solutions outside the recurrence family.

### 4.4 Canonical operator

**Corollary 3 (Canonical K).** All constructions differing in starting vertex,
orientation, or phase seed produce the same K = d₁†⋆₂ d₁.

*Proof.* By Corollary 2, any two constructions differ by d̃₁[f,:] = λ_f d₁[f,:]
with |λ_f| = 1 (since both have unimodular entries). Since
K = Σ_f ⋆₂[f,f] · d₁[f,:]† d₁[f,:], each term acquires |λ_f|² = 1. □

Verified: ‖ΔK‖ < 5.2 × 10⁻¹⁵ under cyclic/reverse face permutations.
Eigenvalues identical (< 10⁻¹²) under 5 random vertex gauge transforms.

### 4.5 Gauge kernel

**Proposition 2 (Kernel dimension).** For generic k (exp(ik_j L_j) ≠ 1 for
all j), dim ker K(k) = |V|.

The exact complex 0 → C⁰ → C¹ → C² → C³ → 0 with operators d₀(k), d₁(k),
d₂(k) is a cochain complex with coefficients in the rank-1 local system L_k
on T³ defined by the representation π₁(T³) → U(1) with monodromy exp(ik_j L_j)
along the j-th cycle. It therefore computes the cellular cohomology
H*(T³, L_k).

*Proof.* Three steps:
1. d₀(k) injective → im d₀ has dimension V.
2. Exactness → ker K ⊇ im d₀, so dim ker K ≥ V.
3. Upper bound: T³ = S¹ × S¹ × S¹, and L_k is a rank-1 local system with
   monodromy λ_j = exp(ik_j L_j) along the j-th cycle. For S¹ with
   monodromy λ ≠ 1, the cohomology H⁰(S¹, L) = ker(λ - 1) = 0 and
   H¹(S¹, L) = coker(λ - 1) = 0, so S¹ is acyclic. By the Künneth formula
   [Bott-Tu, Ch. I], H^p(T³, L_k) = ⊗_j H^*(S¹_j, L_j) = 0 for all p
   when all three monodromies are nontrivial. Hence ker d₁ = im d₀ and
   dim ker K(k) = dim im d₀(k) = V. □

[Verified: rank chain rank(d₀) = ker(d₁) = ker(K) = V on all structures.
Physical modes have gradient fraction < 1.3 × 10⁻¹³.]

### 4.6 Full cochain complex and cohomology

The same recurrence extends to d₂(k), producing the full exact complex:

    0 → C⁰ →^{d₀(k)} C¹ →^{d₁(k)} C² →^{d₂(k)} C³ → 0

with d₂(k)d₁(k) = 0 at machine precision (10⁻¹⁵ on all 4 structures).

**Remark.** The exact complex computes the twisted de Rham cohomology:
- β(Γ) = (1, 3, 3, 1), consistent with H*(T³).
- β(k≠0) = (0, 0, 0, 0), by Künneth with acyclic local coefficients.

Both verified on all 4 test structures.

At k = 0 (mod reciprocal lattice), all monodromies equal 1 and the cohomology
is H*(T³) = (1,3,3,1). At all other k — including TRIM points where some
monodromies equal -1 — the cohomology vanishes: β = (0,0,0,0). This is
because any nontrivial monodromy (λ ≠ 1, even λ = -1) makes the corresponding
S¹ factor acyclic, and an acyclic factor kills the entire tensor product.
Verified at all 4 TRIM points on Kelvin: X (β₁ = 0), M (β₁ = 0), R (β₁ = 0),
BZ_x ≡ Γ (β₁ = 3, n_zero = 98 > V = 96). The cohomology result is thus
stronger than "generic k": it holds at all k ≠ 0 mod reciprocal lattice,
with no exceptions at high-symmetry points.

**Remark (structure of the obstruction).** The discrete model separates into
three layers: topology (incidence), connection (Bloch phases), and metric
(Hodge stars). Exactness depends only on the first two — this is why it
survives arbitrary vertex perturbation (§4.2) but fails under incorrect
phase assignment (§3). The metric enters only in §7 (Voronoi isotropy).

The ratio d₁_std / d₁_exact is a pure U(1) phase per nonzero entry. Its
product around each face boundary — the curvature of the standard connection
relative to the exact one — has a rigid topological structure: the set of
faces with nonzero curvature depends on the k-direction but not on |k|,
the total flux is linear in |k| with integer coefficients that depend on
the mesh structure, and n_curv > n_spur always (local curvature partially
cancels globally). The exact construction is the unique flat connection
on this bundle. Details in Reproducibility (test_07_curvature).

---

## 5. Spectral consequences

### 5.1 Spectral pollution eliminated in all tested cases

| Structure | n_zero exact | n_zero standard | n_spur |
|-----------|-------------|----------------|--------|
| Kelvin N=2 | 96 = V | 90 | 6 |
| C15 N=1 | 136 = V | 127 | 9 |
| WP N=1 | 46 = V | 43 | 3 |
| SC N=3 | 27 = V | 22 | 5 |

Exact: n_zero = V at all k-points, all directions. Universal.
Standard: n_zero < V, direction-dependent (6 on axis, 14 on body diagonal).

### 5.2 Hodge splitting restored

| | Exact | Standard |
|---|-------|---------|
| Mixed modes (out of 192) | 0 (0%) | 140–177 (73–92%) |
| Max gradient overlap (physical) | 1.3 × 10⁻¹³ | — |
| Mean gradient overlap (spurious) | — | 0.88–0.94 |

On standard DEC, the first non-zero eigenvalue is a spurious mode (88–94%
gradient contamination), not acoustic. The acoustic mode cannot be identified
without exactness.

### 5.3 Dispersion convergence

| N | c²_exact | error | c²_standard |
|---|---------|-------|-------------|
| 2 | 0.999679 | 3.2×10⁻⁴ | 1.55 |
| 3 | 0.999857 | 1.4×10⁻⁴ | 2.48 |
| 4 | 0.999920 | 8.0×10⁻⁵ | 1.93 |
| 5 | 0.999949 | 5.1×10⁻⁵ | 1.03 |

Exact: convergence exponent p = 2.00 (O(1/N²) dispersion).
Standard: c² oscillates (mean |c²-1| = 0.75, spread 1.45). Non-convergent.
At N=5 the standard value c²_std = 1.03 happens to be near 1, but the
corresponding eigenmode has 91% gradient contamination — it is a spurious
mode that accidentally coincides with the acoustic frequency, not a correctly
computed acoustic mode.

Convergence confirmed on second mesh family (SC N=3..6): same exponent p=2.00.

[FIGURE 2: Convergence plot — |c²-1| vs N, exact (p=2.00) vs standard
(oscillating). Log-log, Kelvin + SC.]

### 5.4 Universality

Exactness + correct kernel + zero gradient contamination verified on all 5
structure types including 10/10 random Voronoi seeds.

[FIGURE 3: Band structure Γ-X-M-R-Γ for Kelvin N=2. Three panels:
(a) Exact: clean bands, twofold acoustic degeneracy.
(b) Standard: spurious bands collapse toward zero.
(c) ‖d₁d₀‖ along BZ path.]

[TABLE: Universality — ‖d₁d₀‖, n_zero, c², gradient overlap on all 5 structures.]

---

## 6. Consequences of inexactness

Three structural results explaining the spectral pathology of §5. These
connect the curvature structure of §4.6 (local phase inconsistency on each
face) to the global spectral consequences (rank loss, hybridization, eigenvalue
redistribution).

### 6.1 Rank–pollution identity

**Proposition 3 (Rank–pollution identity).**
rank(d₁_std) − rank(d₁_exact) = n_spur.

*Argument.* Both d₁_exact and d₁_std have identical sparsity patterns (nonzero
exactly on (f,e) with e ∈ ∂f) and unimodular entries (|d₁[f,e]| = 1). They
differ only in phases. On the exact complex, ker K = ker d₁ has dimension V
(Proposition 2). On the standard complex, ker K_std has dimension V − n_spur
(n_spur eigenvalues are expelled from the kernel into the physical spectrum).
Since ker K = ker d₁ when ⋆₂ > 0, the kernel dimension change equals the
rank change: dim ker d₁_std − dim ker d₁_exact = −n_spur, hence
rank(d₁_std) − rank(d₁_exact) = n_spur. □

Verified on all structures and all k-directions:

    Kelvin [100]: 102 - 96 = 6 = n_spur
    Kelvin [110]: 108 - 96 = 12 = n_spur
    Kelvin [111]: 110 - 96 = 14 = n_spur

### 6.2 Hodge hybridization at random baseline

On exact DEC: gauge modes have gradient overlap 1.000, physical modes 0.000.
On standard DEC: spurious modes have overlap 0.514 (Kelvin), 0.524 (C15),
0.513 (WP) — within 3% of the random baseline nV/nE = 0.500. This means
the spurious modes show no statistically significant deviation from random
gradient/curl content. The failure is in the deep non-perturbative regime,
not a small correction.

### 6.3 Trace conservation

tr(M⁻¹K) is identical for exact and standard DEC:

    rel_diff < 1.78 × 10⁻¹⁶ across 5 k-directions

This is a topological invariant: |d₁[f,e]| = 1 for both constructions on the
same support, so the trace depends only on the incidence structure and Hodge
stars, not on the phases. Consequence: spurious eigenvalues cannot simply
appear from nothing — they must be compensated by shifts in physical
eigenvalues. The spectral damage is redistributed, not created.

---

## 7. Application: Voronoi tessellations

The exact construction enables identification of geometries where the
curl–curl operator has optimal spectral properties.

### 7.1 Metric isotropy

**Structural property (Voronoi isotropy).** On periodic Voronoi tessellations
of T³, we observe that the discrete metric tensors satisfy G = Vol · I and
H = Vol · I. The identity can be understood via the divergence theorem on
dual cells (for G) and primal cells (for H), exploiting the perpendicularity
of Voronoi edges and dual faces. A complete proof is provided in the
Supplementary Material (Section S1).

**Verified:** ‖G/Vol - I‖ < 10⁻¹⁵ on cubic structures, < 5 × 10⁻¹² on
random Voronoi (5 seeds). Robust under near-degeneracy: with generators
clustered within δ = 0.01 (min ⋆₁ = 3.5 × 10⁻⁵), the identity holds at
3.5 × 10⁻¹².

### 7.2 Exact dispersion

Metric isotropy combined with Theorem 1 (exactness) yields:

    ω² = |k|² + O(|k|⁴)

The argument proceeds by Schur complement reduction of K(k) onto the 3D
harmonic subspace at Γ; exactness ensures this subspace is well-defined
(§4.5), and G = H = Vol·I makes the effective operator proportional to |k|²I.
Full derivation in [Paper A2 / SM].

**Verified:** c² = 0.9993–0.9997 on cubic structures, 0.9995–0.9996 on
random Voronoi (5 seeds). Standard DEC: c²_std = 0.52–1.55 (spurious, not
acoustic). Both metric tensors matter: perturbing ⋆₂ gives Δc² = 0.024,
perturbing ⋆₁ gives Δc² = 0.012. Neither alone is sufficient.

Full spectral analysis, including the converse and moment conservation, in [A2].

---

## 8. Discussion

### Comparison with finite elements

A quantitative comparison with Nédélec edge elements on the same eigenvalue
problems is not included in this work. Such a comparison — measuring L² and
H(curl) errors on standard benchmark problems — would situate the DEC
construction relative to established FEEC methods and is planned for a
companion study.

### Limitations

The Hodge stars are diagonal (circumcentric dual); non-diagonal stars would
require modification of the assembly but not the exactness construction. The
mesh must have valid polyhedral structure; for random Voronoi, sufficient point
density is required (≥50 points per periodic box gives >80% valid tessellations).
The construction is lowest-order (piecewise-constant cochains); extension to
higher-order DEC would require a corresponding generalization of the recurrence.
The observed O(1/N²) convergence rate is consistent with a lowest-order
discretization but a formal error estimate is not attempted.

### Open questions

- Non-flat connections (magnetic flux): the recurrence requires H_f = 1 (flat
  holonomy). For F ≠ 0, a modified recurrence is an open mathematical problem.
- Non-abelian gauge: matrix-valued d₀/d₁ infrastructure not yet available.
- Non-T³ topology: lens spaces, Seifert manifolds require different mesh builders.
- Converse of Voronoi optimality (c²=1 implies G=Vol·I) and moment conservation
  → [Paper A2].
- 2D validation with photonic crystal benchmarks → [Paper B].

## 9. Conclusion

The standard Bloch-periodic extension of DEC does not preserve the discrete
exact sequence on unstructured polyhedral meshes. This structural failure
produces spectral pollution, destroys the Hodge decomposition, and prevents
convergence. A face-boundary recurrence restores exactness by deriving curl
phases from those in the gradient, producing a canonical curl–curl operator
with correct gauge kernel. The construction extends to the full cochain complex
(d₂d₁ = 0) and computes the correct twisted de Rham cohomology. On periodic
Voronoi tessellations, exactness combined with metric isotropy yields exact
isotropic wave speed. The failure is algebraic — independent of mesh quality,
resolution, or symmetry — and is confirmed on five qualitatively different
polyhedral complexes including random Voronoi meshes with no symmetry.

---

## Acknowledgments

[TODO]

## References

[TODO — target 25-35 references. Key citations:]

### FEEC
- Arnold, Falk, Winther (2006). Acta Numerica.
- Arnold, Falk, Winther (2010). Bull. AMS.
- Christiansen (2008). Numerische Mathematik.

### DEC
- Hirani (2003). PhD thesis, Caltech.
- Desbrun, Hirani, Leok, Marsden (2005). arXiv.

### DEC for Maxwell
- Schulz et al. (2018). J. Comp. Phys.
- Mönkölä (2023). Int. J. Numer. Meth. Eng.
- Bossavit (1998). Computational Electromagnetism. Academic Press.

### Spectral pollution
- Boffi (2010). Acta Numerica.
- Boffi, Conforti, Gastaldi (2006). Math. Comp.

### Related methods
- Di Pietro, Droniou (2020). Springer.
- Hiptmair (2002). Acta Numerica.
- Teixeira (1999). IEEE Trans.

### Topology
- Bott, Tu (1982). Differential Forms in Algebraic Topology. Springer.

---

## Reproducibility

All computations use Python 3.9+ with NumPy 1.24+ and SciPy 1.11+. Source code
and complete test suite available at [repository URL].

**Test suite.** 57 tests organized in 7 files mirroring the paper structure:

| File | §  | Tests | What |
|------|-----|-------|------|
| test_01_setup | 2 | 6 | Structure validation, Euler χ, Hodge stars |
| test_02_standard_fails | 3 | 12 | Exactness failure, Cor 1, scaling, sensitivity |
| test_03_recurrence | 4 | 14 | Theorem 1, holonomy, uniqueness, kernel, cohomology, iff |
| test_04_spectral_gains | 5 | 8 | Pollution, Hodge splitting, convergence, universality |
| test_05_inexactness | 6 | 3 | Rank identity, hybridization, trace conservation |
| test_06_voronoi | 7 | 8 | Metric isotropy, c²=1, necessity, near-degeneracy |
| test_07_curvature | 4/6 | 5 | Curvature structure, flux, cancellation |

Total runtime: ~70s on a single core. Each test file includes expected output
in its header docstring for comparison. The mapping from test ID to paper claim
is documented in `tests/tests_map.md`.

**Running.** From the project root:

    python3 run_tests.py                    # all tests
    python3 run_tests.py recurrence voronoi # specific groups

**Structures.** Five test structures are used throughout: SC cubic (N=3),
three periodic Voronoi tessellations (Kelvin/BCC N=2, C15/Laves N=1,
Weaire-Phelan/A15 N=1), and random periodic Voronoi (10 seeds, n=50 cells).
Structure builders are deterministic given seed and size parameters.
