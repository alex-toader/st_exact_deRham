# Exactness-preserving discrete de Rham complexes under Bloch-periodic boundary conditions

**Alexandru Toader**
Independent researcher, Buzău, Romania
toader_alexandru@yahoo.com

---

## Abstract

The standard Bloch-periodic extension of discrete exterior calculus (DEC)
breaks the exactness of the discrete de Rham sequence on unstructured
polyhedral meshes: d₁(k)d₀(k) ≠ 0 for generic wave vector k. We prove
that this failure is structural — no per-edge phase assignment can restore
it — and construct an explicit face-boundary recurrence that produces a
discrete curl d₁(k) satisfying d₁(k)d₀(k) = 0 for all k. The construction
is unique up to a scalar per face and defines a canonical curl–curl operator
whose gauge kernel has dimension |V|. This can be interpreted as constructing
the unique flat discrete connection on the rank-1 local system defined by
Bloch phases. The gauge kernel dimension is determined via the Künneth
formula on T³ with local coefficients, and the full twisted cochain complex
computes the correct de Rham cohomology. Numerical verification on five
polyhedral complexes confirms that exactness eliminates spectral pollution
and restores the Hodge splitting.

**2020 Mathematics Subject Classification:** 65N30, 58A14, 65N25.

**Keywords:** discrete exterior calculus, de Rham complex, Bloch-periodic
boundary conditions, spectral pollution, flat connection.

---

## 1. Introduction

The discrete exterior calculus (DEC) provides a coordinate-free
discretization of differential forms on cell complexes [5, 9]. For
Maxwell's equations on periodic structures, the electric field is
represented as a 1-cochain and the magnetic flux as a 2-cochain, with
discrete gradient d₀ and curl d₁ acting between cochain spaces. The
curl–curl eigenvalue problem

    K a = ω² M a,    K = d₁† ⋆₂ d₁,    M = ⋆₁,

yields electromagnetic eigenfrequencies. Boffi [3] established that
exactness of the discrete complex — the property d₁d₀ = 0 — is both
necessary and sufficient for spurious-free eigenvalue computations. When
exactness holds, the gauge kernel ker K coincides with im d₀, and the
positive-definite Hodge star ⋆₂ ensures discrete compactness of the
physical subspace.

On structured grids (Yee lattice), exactness under Bloch-periodic
boundary conditions holds automatically by tensor-product topology
[10, 11]. Existing DEC-based Maxwell computations on periodic domains
use structured meshes [10, 11]; the unstructured polyhedral case has
not been addressed. The finite element
exterior calculus (FEEC) framework [1, 2, 8] provides the theoretical
foundation for structure-preserving discretizations on simplicial meshes,
but does not explicitly address the phase consistency required by
nontrivial Bloch parameters on polyhedral DEC. The discrete de Rham
(DDR) method [6] achieves exactness on polyhedral meshes but has not
been formulated with Bloch periodicity. To our knowledge, no construction
has been reported that produces an exact discrete de Rham complex under
Bloch-periodic boundary conditions on unstructured polyhedral meshes.

The core observation is that Bloch boundary conditions define a rank-1
local system L_k on the flat torus T³, parametrized by the wave vector
k ∈ R³. The discrete operators must form a cochain complex with
coefficients in L_k. The standard DEC approach assigns Bloch phases
independently per edge — this is a discrete connection with nonzero
curvature. Our construction derives all curl phases from the gradient
via a face-boundary recurrence, producing the unique flat connection
that preserves exactness.

We establish four results:
(i) the standard Bloch extension breaks exactness on any unstructured
mesh, and no per-edge phase fix exists (Proposition 1, Corollary 1);
(ii) an explicit recurrence constructs d₁(k) with d₁(k)d₀(k) = 0 for
all k (Theorem 1);
(iii) the construction is unique up to per-face gauge and defines a
canonical curl–curl operator (Theorem 1);
(iv) the gauge kernel has dimension |V| for generic k, computed via the
Künneth formula on T³ with coefficients in L_k (Proposition 2).

The construction is verified on five periodic polyhedral complexes
including random Voronoi meshes with no symmetry (§5).¹

¹ Source code and test suite: github.com/alex-toader/st_exact_deRham.

---

## 2. Periodic cell complexes and Bloch operators

### 2.1. Cell complex on the torus

Let Λ = L₁Z × L₂Z × L₃Z be a lattice in R³ with fundamental domain
Ω = [0,L₁) × [0,L₂) × [0,L₃). A periodic cell complex (V, E, F, C)
on the flat torus T³ = R³/Λ consists of vertices V, oriented edges E,
oriented faces F, and 3-cells C. Edges connecting vertices across the
periodic boundary carry a lattice shift vector nₑ ∈ Z³; interior edges
have nₑ = 0.

The topological incidence matrices d₀ ∈ Z^{|E|×|V|} and d₁ ∈ Z^{|F|×|E|}
encode boundary relations and satisfy d₁d₀ = 0.

### 2.2. Bloch-twisted operators

For wave vector k ∈ R³, the Bloch-twisted gradient is

    (d₀(k) u)ₑ = e^{ik·nₑL} · u_{head(e)} − u_{tail(e)}.

The standard Bloch curl applies per-edge phases independently:

    d₁_std(k)[f, e] = d₁[f, e] · e^{ik·nₑL}.

At k = 0, both reduce to the topological operators.

### 2.3. Hodge stars and curl–curl operator

The diagonal Hodge stars ⋆₁[e,e] = |σ*_e|/|σ_e| and ⋆₂[f,f] = |σ*_f|/|σ_f|
encode metric information from the circumcentric dual. The curl–curl
eigenvalue problem is K a = ω² ⋆₁ a with K = d₁† ⋆₂ d₁.

### 2.4. Test structures

The construction is verified on five periodic polyhedral complexes: the
SC cubic lattice (N=3), three periodic Voronoi tessellations (Kelvin/BCC
N=2, C15/Laves N=1, Weaire-Phelan/A15 N=1), and random Voronoi meshes
(50 cells, 10 seeds). The Voronoi tessellations have 4-valent vertices
and 3 faces per edge (Plateau structure).

---

## 3. Failure of the standard construction

**Proposition 1** (Structural failure).
*Let f be a face whose boundary contains two edges eₐ, e_b sharing a
vertex, with distinct lattice shifts nₐ ≠ n_b. Then*

    d₁_std(k) d₀(k) ≠ 0

*for all k outside the hyperplane {k : e^{ik·(nₐ−n_b)·L} = 1}. The failure
occurs on an open dense subset of R³.*

*Proof.* At the shared vertex v, the two contributions to (d₁d₀)[f,v]
carry phases e^{ik·nₐL} and e^{ik·n_bL}. Their cancellation requires
e^{ik·(nₐ−n_b)L} = 1, which defines a hyperplane in k-space. This holds
at k = 0 but fails for generic k when nₐ ≠ n_b. □

**Corollary 1** (No per-edge fix).
*Let φ: E × R³ → C* assign a phase to each edge, independent of the face
context. Then d̃₁(k)d₀(k) ≠ 0 for generic k on any mesh with faces whose
boundary edges carry distinct lattice shifts.*

*Proof.* A per-edge multiplicator d̃₁[f,e] = d₁[f,e] · φ(e,k) assigns
the same phase to edge e in every face containing it. But the exactness
condition d₁d₀ = 0 at face f determines the ratio φᵢ/φ₀ = e^{ik·SᵢL}
where Sᵢ ∈ Z³ is the net lattice shift along the partial boundary path
from e₀ to eᵢ (as shown by the recurrence constructed in §4). Since
Sᵢ depends on the position of edge eᵢ within face f — not on eᵢ alone —
two faces sharing an edge e at different boundary positions require
different phases for e. A per-edge function cannot produce this
face-dependent variation. □

On all five test structures, ‖d₁_std(k) d₀(k)‖ ∈ [5, 12] at generic
k-points; the exact construction of §4 achieves ‖d₁(k) d₀(k)‖ < 10⁻¹⁵
(see §5 for details).

---

## 4. Exactness-preserving construction

**Theorem 1** (Exactness-preserving Bloch-twisted complex).
*Let (V, E, F, C) be a periodic polyhedral cell complex on T³ and
k ∈ R³. There exists a local operator d₁(k): C¹ → C² with support
on face boundaries satisfying*

    d₁(k) d₀(k) = 0    for all k ∈ R³.

*For each face f with ordered boundary edges e₀, ..., eₙ₋₁ and incidence
signs σᵢ = d₁[f, eᵢ] ∈ {±1}, the operator is given by the recurrence:*

    (i)   φ₀ = 1,
    (ii)  φᵢ = −σᵢ₋₁ φᵢ₋₁ · d₀(k)[eᵢ₋₁, vᵢ] / (σᵢ · d₀(k)[eᵢ, vᵢ]),
    (iii) d₁(k)[f, eᵢ] = σᵢ φᵢ.

*The construction is unique up to a nonzero scalar per face: any local
operator with d₁d₀ = 0 supported on ∂f satisfies d₁[f,:] = λ_f · d₁^{ex}[f,:]
for some λ_f ∈ C* with |λ_f| = 1. The solution space is one-dimensional
per face.*

*The curl–curl operator K(k) = d₁(k)† ⋆₂ d₁(k) is independent of the
choice of starting vertex, boundary orientation, or per-face gauge λ_f.*

*Proof.* The condition d₁(k)d₀(k) = 0 at face f with n boundary edges
decomposes into n vertex equations. At vertex vᵢ for i = 1, ..., n−1:

    d₁[f, eᵢ₋₁] · d₀(k)[eᵢ₋₁, vᵢ] + d₁[f, eᵢ] · d₀(k)[eᵢ, vᵢ] = 0.

Substituting d₁[f, eᵢ] = σᵢ φᵢ gives the recurrence (ii), which
determines φ₁, ..., φₙ₋₁ from φ₀ = 1. The remaining equation at v₀
involves d₁[f, eₙ₋₁] and d₁[f, e₀] and is satisfied if and only if
the product of all recurrence ratios around the face boundary equals 1.
One verifies that this product equals the holonomy H_f = ∏ᵢ e^{ik·nᵢL} around ∂f,
which we now show is trivial.

**Lemma** (Flat holonomy). *For any face f, the holonomy product around
the boundary is trivial:*

    H_f = ∏ᵢ e^{ik·nᵢL} = e^{ik·(Σᵢ nᵢ)L} = 1.

*Proof.* Choose a lift of f to the universal cover R³. The lifted boundary
is a closed path in R³ bounding a contractible polygon. Consecutive edges
concatenate without jumps, so the total lattice translation is Σᵢ nᵢ = 0.
Hence H_f = e^{ik·0} = 1. □

The lemma ensures that the n vertex equations are not independent: they
sum to zero by the topological identity d₁d₀ = 0 at k = 0. Any n−1 of
them imply the remaining one. The recurrence therefore closes consistently
for every face.

For uniqueness: the system is a first-order recurrence with n−1 equations
in n unknowns (the entries d₁[f, eᵢ]), hence the solution space is
one-dimensional per face.

For K canonical: two constructions differing by d̃₁[f,:] = λ_f d₁[f,:]
with |λ_f| = 1 contribute |λ_f|² = 1 to each term of
K = Σ_f ⋆₂[f,f] · d₁[f,:]† d₁[f,:]. □

**Remark** (Geometric viewpoint). Bloch boundary conditions define a
rank-1 local system L_k on T³ with monodromy e^{ikⱼLⱼ} along the j-th
cycle. The standard construction is a discrete connection on L_k with
nonzero curvature — the phase inconsistency of Proposition 1 is
precisely this curvature. The recurrence of Theorem 1 produces the unique
flat connection on L_k induced by d₀(k). Exactness of the twisted
cochain complex is equivalent to flatness of this discrete connection.

**Proposition 2** (Kernel dimension).
*For generic k (i.e., e^{ikⱼLⱼ} ≠ 1 for all j = 1,2,3),*

    dim ker K(k) = |V|    and    ker K(k) = im d₀(k).

*Proof.* Three steps.

(1) d₀(k) is injective: if d₀(k)u = 0, then for each edge e with
lattice shift nₑ, we have e^{ik·nₑL} u_{head} = u_{tail}. Iterating
along any cycle generating π₁(T³) with monodromy λⱼ = e^{ikⱼLⱼ} ≠ 1
forces u = 0. Hence dim im d₀(k) = |V|.

(2) Exactness (Theorem 1) gives im d₀ ⊆ ker d₁. Since ⋆₂ > 0,
ker K = ker d₁, so ker K ⊇ im d₀ and dim ker K ≥ |V|.

(3) For the reverse inequality, T³ = S¹ × S¹ × S¹ and L_k restricts
to a rank-1 local system on each factor with monodromy λⱼ ≠ 1. On S¹
with nontrivial monodromy, H⁰ = ker(λ−1) = 0 and H¹ = coker(λ−1) = 0,
so each factor is acyclic. By the Künneth formula [4, 7], H^p(T³, L_k) = 0
for all p, hence ker d₁ = im d₀ and dim ker K = |V|. □

The same recurrence extends to d₂(k): C² → C³, producing d₂(k)d₁(k) = 0
at machine precision. The full complex computes the correct twisted
de Rham cohomology: β(Γ) = (1, 3, 3, 1) consistent with H*(T³), and
β(k) = (0, 0, 0, 0) for k ≠ 0 mod reciprocal lattice.

---

## 5. Numerical verification

Table 1 compares the gauge kernel dimension at k = (0.1, 0, 0)·2π/L:

| Structure | |V| | n_zero (exact) | n_zero (standard) | n_spur |
|-----------|-----|---------------|------------------|--------|
| Kelvin N=2 | 96 | 96 | 90 | 6 |
| C15 N=1 | 136 | 136 | 127 | 9 |
| WP N=1 | 46 | 46 | 43 | 3 |
| SC N=3 | 27 | 27 | 22 | 5 |

On the exact complex, the gauge kernel has dimension |V| at every tested
k-point and direction, confirming Proposition 2. On the standard complex,
n_spur eigenvalues are expelled from the kernel into the physical spectrum,
producing spurious modes.

The Hodge splitting is fully restored by the exact construction: out of
192 eigenmodes on Kelvin N=2, zero show gradient contamination (gradient
overlap < 10⁻¹³). On the standard complex, 73–92% of modes are hybridized,
with gradient overlap 0.88–0.94 — consistent with the random baseline
|V|/|E| = 0.5.

The phase velocity converges at rate O(1/N²) on Kelvin (N = 2, ..., 5)
with exponent p = 2.00, confirmed independently on SC (N = 3, ..., 6).
The standard construction does not converge: c² oscillates between 1.03
and 2.48.

Figure 1 shows the band structure along Γ–X–M–R–Γ for Kelvin N=2.
Panel (a): the exact construction produces clean bands with twofold
acoustic degeneracy. Panel (b): the standard construction produces
spurious bands collapsing toward zero frequency.

[FIGURE 1: Band structure Γ–X–M–R–Γ, Kelvin N=2.
(a) Exact DEC. (b) Standard DEC.]

Exactness and correct kernel dimension are verified on all five structure
types, including 10 out of 10 random Voronoi seeds.

---

## 6. Application to Voronoi tessellations

**Remark** (Voronoi optimality). On any periodic Voronoi tessellation of T³,
the discrete metric tensors satisfy G = Vol · I and H = Vol · I (where Vol denotes the volume of
the fundamental domain Ω), with
G_{ij} = Σₑ (⋆₁)ₑ (ℓₑ)ᵢ (ℓₑ)ⱼ and H_{ij} = Σ_f (⋆₂)_f (Aᶠ)ᵢ (Aᶠ)ⱼ,
where ℓₑ is the edge direction vector (unit tangent × length) and Aᶠ
is the face area vector (unit normal × area), summed over
boundary-crossing edges and faces, respectively. This
follows from the divergence theorem applied to dual cells (for G) and
primal cells (for H), exploiting the perpendicularity of Voronoi edges
to dual faces (details in [12]). Verified: ‖G/Vol − I‖ < 10⁻¹⁵ on cubic
structures and < 5 × 10⁻¹² on random Voronoi (5 seeds).

Combined with Theorem 1, metric isotropy yields the exact acoustic
dispersion relation ω² = |k|² + O(|k|⁴) via Schur complement reduction
of K(k) onto the three-dimensional harmonic subspace at Γ. The effective
operator is proportional to |k|²I when G = H = Vol · I. Verified:
c² ∈ [0.9993, 0.9997] on cubic structures and [0.9995, 0.9996] on random
Voronoi seeds. On the standard complex, c² ranges from 0.52 to 1.55.

---

## 7. Conclusion

The standard Bloch-periodic extension of DEC does not preserve the exact
sequence on unstructured polyhedral meshes — a structural failure, not a
numerical one. The face-boundary recurrence of Theorem 1 produces a
canonical exact complex, unique up to per-face gauge, with kernel
dimension |V| and correct twisted de Rham cohomology. Extensions to
non-flat connections (nonzero magnetic flux), non-abelian gauge groups,
and non-T³ topologies remain open.

---

## References

[1] D.N. Arnold, R.S. Falk, R. Winther, *Finite element exterior calculus, homological techniques, and applications*, Acta Numer. **15** (2006), 1–155.

[2] D.N. Arnold, R.S. Falk, R. Winther, *Finite element exterior calculus: from Hodge theory to numerical stability*, Bull. Amer. Math. Soc. **47** (2010), 281–354.

[3] D. Boffi, *Finite element approximation of eigenvalue problems*, Acta Numer. **19** (2010), 1–120.

[4] R. Bott, L.W. Tu, *Differential Forms in Algebraic Topology*, Springer, 1982.

[5] M. Desbrun, A.N. Hirani, M. Leok, J.E. Marsden, *Discrete exterior calculus*, arXiv:math/0508341, 2005.

[6] D.A. Di Pietro, J. Droniou, *The Hybrid High-Order Method for Polytopal Meshes*, Springer, 2020.

[7] A. Hatcher, *Algebraic Topology*, Cambridge University Press, 2002.

[8] R. Hiptmair, *Finite elements in computational electromagnetism*, Acta Numer. **11** (2002), 237–339.

[9] A.N. Hirani, *Discrete Exterior Calculus*, PhD thesis, Caltech, 2003.

[10] S. Mönkölä, J. Räbinä, T. Rossi, *Discrete exterior calculus for photonic crystal waveguides*, Int. J. Numer. Methods Eng. **124** (2023), 1035–1054.

[11] F.L. Teixeira, W.C. Chew, *Lattice electromagnetic theory from a topological viewpoint*, J. Math. Phys. **40** (1999), 169–187.

[12] A. Toader, *Metric isotropy and spectral optimality on periodic Voronoi tessellations*, in preparation.
