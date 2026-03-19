# Abstract drafts

## Draft 1 (17 Mar 2026)

The discrete exterior calculus (DEC) on periodic polyhedral complexes is widely
used for eigenvalue problems in computational electromagnetics. Under
Bloch-periodic boundary conditions, however, the standard construction of the
discrete curl operator fails to preserve the de Rham exact sequence:
$d_1(\mathbf{k})\,d_0(\mathbf{k}) \neq 0$ for generic wave vector $\mathbf{k}$
on any unstructured mesh with faces whose boundary edges carry distinct lattice
shifts. We prove that no per-edge phase assignment can restore exactness, and
construct an alternative $d_1(\mathbf{k})$ via a face-boundary recurrence that
enforces $d_1(\mathbf{k})\,d_0(\mathbf{k}) = 0$ for all $\mathbf{k}$. The
construction is explicit, costs $O(\sum n_f)$, and is unique up to a $U(1)$
phase per face; the resulting curl--curl operator is canonical. We show that
the gauge kernel has dimension $|V|$ by computing the twisted de Rham
cohomology $H^p(T^3, \mathcal{L}_\mathbf{k})$ via the Künneth formula.
Without exactness, we identify three structural consequences: the pollution
count equals the rank excess of $d_1$ (an algebraic identity), the Hodge
splitting degrades to the random-baseline overlap $|V|/|E|$, and the total
eigenvalue trace is conserved between shifted physical modes and spurious
modes. Numerical experiments on five polyhedral complexes---including Voronoi
foams and random tessellations---confirm that exactness eliminates spectral
pollution, restores the discrete Hodge decomposition, and yields second-order
dispersion convergence. On periodic Voronoi tessellations, the metric tensors
satisfy $G = \mathrm{Vol}\cdot I$ and the construction gives exact isotropic
wave speed $\omega^2 = |\mathbf{k}|^2 + O(|\mathbf{k}|^4)$.

---

Word count: ~220. Target: 150-200 for SIAM.

## Notes on Draft 1

- ~220 words, target 150-180
- Too dense: Künneth, all 3 W3 results, Voronoi formula
- Reviewer feedback: cut topology from abstract, reduce §6 to one sentence

---

## Draft 2 (17 Mar 2026) — superseded by Draft 3

[kept for reference]

---

## Draft 3 (17 Mar 2026)

Discrete exterior calculus (DEC) is often described as free of spurious
modes in the curl--curl eigenvalue problem. We show that this property fails
under Bloch-periodic boundary conditions: the standard discrete curl satisfies
$d_1(\mathbf{k})\,d_0(\mathbf{k}) \neq 0$ for generic wave vector
$\mathbf{k}$ on any unstructured polyhedral mesh with faces whose boundary
edges carry distinct lattice shifts. No per-edge phase assignment can restore
exactness. We construct an alternative $d_1(\mathbf{k})$ via a face-boundary
recurrence that enforces $d_1(\mathbf{k})\,d_0(\mathbf{k}) = 0$ for all
$\mathbf{k}$. The construction is explicit, unique up to per-face gauge,
and produces a canonical curl--curl operator with gauge kernel of dimension
$|V|$. Without exactness, the pollution count equals the rank excess of
$d_1$ and the Hodge decomposition degrades to random-baseline hybridization.
Experiments on five polyhedral complexes confirm that exactness eliminates
spectral pollution, restores the Hodge splitting, and yields second-order
dispersion convergence. On periodic Voronoi tessellations, exactness combined
with metric isotropy gives $\omega^2 = |\mathbf{k}|^2 + O(|\mathbf{k}|^4)$.

## Notes on Draft 3

- ~175 words ✓
- NEW: Opens with contrast — "often described as free of spurious modes... we show this fails"
- Künneth out (stays in §4)
- §6 = one sentence (rank excess + hybridization)
- Voronoi = last sentence
- Reviewer immediately sees: (1) what community believes, (2) why it's wrong, (3) the fix
