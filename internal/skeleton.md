# SKELETON v3 — CR Math Note (FINAL)
# "Exactness-preserving discrete de Rham complexes under Bloch-periodic boundary conditions"
# Target: Comptes Rendus Mathématique, ~8 pages, ~50k caractere

---

## §1. Introduction [~1 pag]

**Paragraf 1 — Context:**
DEC pe complexe poliedrale: d₀, d₁, secvența d₁d₀=0.
Maxwell eigenvalue problem → exactitatea = necesară și suficientă (Boffi 2010).

**Paragraf 2 — Gap:**
Sub Bloch BC pe mesh nestructurate: standard sparge exactitatea.
FEEC [Arnold-Falk-Winther] = simplicial. DEC existent = Yee structurat.
Nicio construcție exactă pe mesh poliedrale general sub Bloch.

**Paragraf 3 — Viewpoint + Contribuție:**
Bloch BC = sistem local rang-1 L_k pe T³. Construcția standard = conexiune
cu curbură nenulă. Noi: unica conexiune plată care păstrează exactitatea.

Contribuții:
(i)   eșec structural al standardului (Prop. 1)
(ii)  construcție explicită prin recurență (Thm. 1)
(iii) unicitate + K canonic (Thm. 1)
(iv)  dim ker K = |V| via Künneth (Prop. 2)

Footnote: Code available at github.com/alex-toader/st_exact_deRham.

---

## §2. Setup [~1 pag]

**§2.1 Complex celular pe T³.**
(V, E, F, C), orientări, shift-uri nₑ ∈ Z³.
Incidență: d₀ ∈ R^{|E|×|V|}, d₁ ∈ R^{|F|×|E|}. Identitate d₁d₀ = 0 la k=0.

**§2.2 Operatori Bloch.**
d₀(k): (d₀(k)u)ₑ = e^{ik·nₑL} u_{head} − u_{tail}. Formulă explicită.
d₁_std(k): d₁_std[f,e] = d₁[f,e] · e^{ik·nₑL}. Faze independente per-muchie.

**§2.3 Hodge stars.**
⋆₁, ⋆₂ diagonale (circumcentric). K = d₁†⋆₂d₁, M = ⋆₁.

**§2.4 Structuri test (o frază).**
"Verified on five periodic polyhedral complexes: SC cubic lattice, three
Voronoi tessellations (Kelvin/BCC, C15/Laves, Weaire-Phelan/A15), and
random Voronoi meshes (50 cells, 10 seeds)."

---

## §3. Eșecul construcției standard [~1 pag]

**Propoziție 1 (Structural failure).**
If face f has ≥2 boundary edges with distinct lattice shifts nₐ ≠ n_b
sharing a vertex, then d₁_std(k)d₀(k) ≠ 0 for all k outside the hyperplane
{e^{ik·(nₐ−n_b)L} = 1}. The failure occurs on an open dense set of k.

*Proof (4-5 rânduri):*
La vârful comun v, contribuțiile la (d₁d₀)[f,v] poartă faze e^{ik·nₐL}
și e^{ik·n_bL}. Anularea cere e^{ik·(nₐ−n_b)L}=1 — hiperplan.
Eșuează generic când nₐ ≠ n_b. □

**Corolar 1 (No per-edge fix).**
No function φ: E×R³ → C* can make d̃₁(k)d₀(k) = 0 generically.
The obstruction is per-face (2-cochain), not per-edge (1-cochain).

*Proof (3-4 rânduri):*
By uniqueness (Thm. 1), exact phases are face-dependent. A per-edge
multiplicator assigns the same factor regardless of face. Contradiction. □

Notă numerică (1 frază): "‖d₁_std d₀‖ ∈ [5, 12] on all tested structures;
the exact construction achieves 10⁻¹⁵ (§5)."

---

## §4. Construcția exactă [~2.5 pag] — NUCLEUL

**Teorema 1 (Exactness-preserving Bloch-twisted complex).**

*Enunț:*
Let (V,E,F,C) be a periodic polyhedral cell complex on T³ and k ∈ R³.
There exists a local operator d₁(k): C¹ → C² supported on face boundaries
satisfying d₁(k)d₀(k) = 0 for all k ∈ R³.

Construction: for face f with ordered boundary edges e₀,...,eₙ₋₁:
  (i)   φ₀ = 1
  (ii)  φᵢ = −σᵢ₋₁ φᵢ₋₁ · d₀(k)[eᵢ₋₁,vᵢ] / (σᵢ · d₀(k)[eᵢ,vᵢ])
  (iii) d₁(k)[f,eᵢ] = σᵢ φᵢ

Uniqueness: any local d₁ with d₁d₀=0 on ∂f satisfies d₁[f,:] = λ_f · d₁^{ex}[f,:]
with λ_f ∈ C*, |λ_f|=1. The solution space is 1-dimensional per face.

Canonical K: K(k) = d₁(k)†⋆₂d₁(k) is independent of seed, orientation,
or per-face gauge.

*Proof:*

→ d₁d₀=0 la face f dă n ecuații la vârfuri. La vᵢ (i=1,...,n-1):
  d₁[f,eᵢ₋₁]·d₀(k)[eᵢ₋₁,vᵢ] + d₁[f,eᵢ]·d₀(k)[eᵢ,vᵢ] = 0.
  Recurența (ii) rezolvă φ₁,...,φₙ₋₁ din φ₀=1.

→ Ecuația rămasă la v₀ e satisfăcută ⟺ holonomia pe față = 1.

  **Lemma (Flat holonomy).** H_f = ∏ e^{ik·nᵢL} = e^{ik·(Σnᵢ)L} = 1,
  deoarece Σnᵢ = 0 (cale închisă în R³ pe lift contractibil al lui f). □_Lemma

→ Unicitate: recurență de ordinul 1, n-1 ecuații în n necunoscute
  → spațiu 1-dimensional per față.

→ K canonic: |λ_f|²=1 se anulează în d₁†⋆₂d₁. □

**Remark (Geometric viewpoint).**
Bloch BC define a rank-1 local system L_k on T³. The standard construction
is a discrete connection on L_k with nonzero curvature. The recurrence of
Theorem 1 produces the unique flat connection induced by d₀(k). Exactness
of the cochain complex is equivalent to flatness of the connection.

**Propoziție 2 (Kernel dimension).**

*Enunț:*
For generic k (e^{ikⱼLⱼ} ≠ 1 for all j), dim ker K(k) = |V|
and ker K(k) = im d₀(k).

*Proof (3 pași):*
1. d₀(k) injective: d₀(k)u=0 forces uᵢ=0 for all i (iterate on toric
   cycles with monodromies ≠ 1). → dim im d₀(k) = |V|.
2. Exactness (Thm. 1) → im d₀ ⊆ ker d₁ → ker K ⊇ im d₀.
3. Künneth on T³ = S¹×S¹×S¹ with coefficients in L_k: each S¹ factor
   with monodromy λⱼ ≠ 1 is acyclic → H*(T³, L_k) = 0 → ker d₁ = im d₀. □

**Complexul complet (1 frază):**
The same recurrence extends to d₂(k), giving d₂d₁=0. Betti numbers:
β(Γ)=(1,3,3,1), β(k≠0)=(0,0,0,0).

---

## §5. Verificare numerică [~0.75 pag]

**Tabel 1** (compact, 4×3):

| Structure | n_zero exact | n_zero std | n_spur |
|-----------|-------------|-----------|--------|
| Kelvin N=2 | 96 = |V| | 90 | 6 |
| C15 N=1 | 136 = |V| | 127 | 9 |
| WP N=1 | 46 = |V| | 43 | 3 |
| SC N=3 | 27 = |V| | 22 | 5 |

**Hodge splitting (1 frază):**
"On the exact complex, 0 out of 192 eigenmodes show gradient contamination;
on the standard complex, 73-92% of modes are hybridized (gradient overlap
0.88-0.94, consistent with the random baseline |V|/|E|)."

**Convergență (1 frază):**
"Phase velocity converges at rate p=2.00 (Kelvin N=2..5, confirmed on SC);
the standard construction does not converge."

**Figură 1:** Band structure Γ-X-M-R-Γ, Kelvin N=2.
  (a) Exact DEC: clean bands, twofold acoustic degeneracy.
  (b) Standard DEC: spurious bands near zero.

**Universalitate (1 frază):**
"Exactness and correct kernel verified on all five structures including
10/10 random Voronoi seeds."

---

## §6. Aplicație: Voronoi [~0.5 pag]

**Remark (Voronoi optimality).**
On any periodic Voronoi tessellation of T³, the discrete metric tensors
satisfy G = H = Vol·I (via the divergence theorem on dual/primal cells;
details in [A2]). Combined with Theorem 1, this yields
ω² = |k|² + O(|k|⁴) by Schur complement reduction onto the 3D harmonic
subspace at Γ. Verified: c² ∈ [0.9993, 0.9997] on cubic structures and
random Voronoi seeds.

---

## §7. Concluzii [~2-3 fraze]

The standard Bloch-periodic extension of DEC does not preserve the exact
sequence on unstructured polyhedral meshes — a structural, not numerical,
failure. The face-boundary recurrence of Theorem 1 produces a canonical
exact complex, unique up to per-face gauge, with kernel dimension |V| and
correct twisted de Rham cohomology. Extensions to non-flat connections,
non-abelian gauge, and non-T³ topologies remain open.

---

## Referințe [~15 intrări]

[1-2]   Arnold, Falk, Winther (FEEC, 2006 + 2010)
[3]     Boffi (eigenvalue problems, 2010)
[4]     Bossavit (computational EM, 1998)
[5-6]   Bott-Tu / Hatcher (Künneth, algebraic topology)
[7]     Christiansen (stability DEC, 2008)
[8]     Desbrun, Hirani, Leok, Marsden (DEC, 2005)
[9]     Di Pietro, Droniou (DDR, 2020)
[10]    Hiptmair (Maxwell FE, 2002)
[11]    Hirani (DEC thesis, 2003)
[12-14] Mönkölä / Schulz / Teixeira (DEC periodic)
[15]    [A2] Toader, in preparation (Voronoi isotropy + converse)

---

## DECIZII DE DESIGN

### Tăiat față de ESAIM
- ✗ Reproducibility section → footnote GitHub în §1
- ✗ Tabelul cu 57 teste → GitHub only
- ✗ §6 vechi (rank-pollution, hybridization, trace conservation) → eliminat complet
- ✗ §7 vechi (Voronoi lung, sensitivitate metrică) → Remark 0.5 pag
- ✗ Discussion (comparație FE, limitări) → eliminat
- ✗ Supplementary material → referință [A2]
- ✗ Subfigura (c) cu ‖d₁d₀‖ → eliminată
- ✗ Tabelul structuri din §2 → frază în proză

### Ce am luat din review-uri
- ✅ Viewpoint "conexiune plată pe L_k" în intro
- ✅ Remark geometric explicit după Teoremă
- ✅ Teorema: bloc complet (existență+construcție+unicitate+K canonic)
- ✅ Main result first: Teorema înaintea Lemei
- ✅ §3 verificare → trimite la §5
- ✅ Voronoi ca Remark
- ✅ Referințe ~15 (nu 35)
- ✅ Subfigura (c) scoasă, doar (a)+(b)

### Ce NU am luat
- ❌ "Taie §5 la o singură afirmație" — tabel + figură + fraze condensate rămân
- ❌ Varianta 2 abstractă de teoremă — prea opacă
- ❌ "Scoate procentele" — rămân într-o frază

## Nucleul matematic (7 enunțuri)

1. **Propoziție 1** §3 — eșecul structural
2. **Corolar 1** §3 — per-edge fix imposibil
3. **Teorema 1** §4 — recurența (existență + construcție + unicitate + K canonic)
4. **Lemma** §4 — holonomie plată (în proof-ul Teoremei)
5. **Remark** §4 — viewpoint: conexiune plată pe L_k
6. **Propoziție 2** §4 — dim ker = |V| via Künneth
7. **Remark** §6 — Voronoi isotropy + ω²=|k|²+O(|k|⁴)
