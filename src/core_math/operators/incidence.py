"""
Incidence Matrices and Hodge Laplacian
======================================

Pure combinatorics - NO physics assumptions.

DEFINITIONS:
    d₀: E × V  "gradient" - oriented edge-vertex incidence
    d₁: F × E  "curl" - oriented face-edge incidence

    L₁ = d₀d₀ᵀ + d₁ᵀd₁  (Hodge Laplacian on edges)

TRACE IDENTITIES (topology-dependent):
    1. d₁d₀ = 0           (exactness - ALWAYS)
    2. Tr(d₀d₀ᵀ) = 2E     (ALWAYS - each edge has 2 endpoints)
    3. Tr(d₁ᵀd₁) = k·E    (k = faces_per_edge from mesh contract)
       - Surface (2-manifold): k = 2
       - Foam (3D Plateau):    k = 3
    4. Tr(L₁) = (2+k)·E   (sum of above)

CYCLE SPACE - IMPORTANT DISTINCTION:

    There are TWO DIFFERENT "H" concepts:

    1. H₁(surface) = surface homology of the 2D manifold
       For Kelvin cell (topological sphere): H₁(surface) = 0

    2. H = ker(d₀ᵀ) = cycle space of the GRAPH (1-skeleton)
       This is what we compute and use.
       dim(H) = E - V + 1 = 36 - 24 + 1 = 13 for Kelvin

    These are DIFFERENT mathematical objects:
    - H₁(surface) counts holes in the surface (none for sphere)
    - H(graph) counts independent cycles in the edge graph

    We work with H(graph), which has dim = 13 even though
    the surface has no 1-cycles.

    DIMENSION FORMULA:
        dim(H) = E - V + 1 = E - rank(d₀)
        This follows from rank-nullity on d₀ᵀ.

REFERENCE: Discrete Exterior Calculus (Desbrun et al., 2005)
"""

import numpy as np
from collections import deque
from typing import List, Tuple, Dict, Any

from ..spec.constants import EPS_CLOSE


def build_d0(vertices: np.ndarray,
             edges: List[Tuple[int, int]]) -> np.ndarray:
    """
    Build gradient operator d₀: C⁰ → C¹.

    DEFINITION:
        d₀[e, v] = -1 if v is the source of edge e
        d₀[e, v] = +1 if v is the target of edge e
        d₀[e, v] = 0 otherwise

    Convention: for edge (i, j) with i < j, i is source, j is target.

    Args:
        vertices: (V, 3) array (only used for V count)
        edges: list of E tuples (i, j)

    Returns:
        d0: (E, V) dense incidence matrix

    PROPERTY:
        Each row has exactly one -1 and one +1.
        Column sum is zero (edges are oriented consistently).
    """
    V = len(vertices)
    E = len(edges)
    d0 = np.zeros((E, V))

    for e_idx, (i, j) in enumerate(edges):
        d0[e_idx, i] = -1  # source
        d0[e_idx, j] = +1  # target

    return d0


def build_d1(vertices: np.ndarray,
             edges: List[Tuple[int, int]],
             faces: List[List[int]]) -> np.ndarray:
    """
    Build curl operator d₁: C¹ → C².

    DEFINITION:
        d₁[f, e] = +1 if edge e appears in face f boundary with positive orientation
        d₁[f, e] = -1 if edge e appears with negative orientation
        d₁[f, e] = 0 if edge e is not on boundary of f

    Args:
        vertices: (V, 3) array
        edges: list of E tuples (i, j)
        faces: list of F faces, each face is vertex list in CCW order

    Returns:
        d1: (F, E) sparse incidence matrix

    PROPERTY:
        Each column has exactly k non-zero entries, where k = faces_per_edge
        from mesh contract (k=2 for surface, k=3 for foam, k=4 for SC solid).
        This is because each edge bounds exactly k faces.

    FAIL-FAST:
        Raises ValueError if any face segment is not in edge list.
    """
    V = len(vertices)
    E = len(edges)
    F = len(faces)

    # Build edge lookup: (i,j) -> (edge_index, sign)
    edge_dict = {}
    for e_idx, (i, j) in enumerate(edges):
        edge_dict[(i, j)] = (e_idx, +1)  # forward direction
        edge_dict[(j, i)] = (e_idx, -1)  # backward direction

    d1 = np.zeros((F, E))

    for f_idx, face in enumerate(faces):
        n = len(face)
        for k in range(n):
            v1 = face[k]
            v2 = face[(k + 1) % n]
            if (v1, v2) not in edge_dict:
                raise ValueError(f"Face {f_idx} uses segment ({v1},{v2}) which is not in edge list. "
                                f"Face vertices: {face}")
            e_idx, sign = edge_dict[(v1, v2)]
            # Use += and detect duplicate edge in same face (invalid cycle)
            if d1[f_idx, e_idx] != 0:
                raise ValueError(f"Face {f_idx} uses edge {e_idx} twice. "
                                f"This indicates an invalid face cycle. Face vertices: {face}")
            d1[f_idx, e_idx] += sign

    return d1


def build_d2(cell_face_incidence: List[List[Tuple[int, int]]],
             n_faces: int,
             verify: bool = True) -> np.ndarray:
    """
    Build divergence operator d₂: C² → C³ (faces to cells).

    DEFINITION (matches spec/constants.py convention):
        d₂[c, f] = +1 if face f is on boundary of cell c with OUTWARD normal
        d₂[c, f] = -1 if face f is on boundary of cell c with INWARD normal
        d₂[c, f] = 0  if face f is not on boundary of cell c

    Args:
        cell_face_incidence: list of lists, cell_face_incidence[c] = [(face_idx, orientation), ...]
                             where orientation = +1 (outward) or -1 (inward)
        n_faces: total number of faces (F)
        verify: if True (default), verify each face appears in exactly 2 cells with sum=0

    Returns:
        d2: (C, F) incidence matrix

    PROPERTY:
        Each face appears in exactly 2 cells (for closed foam), with opposite orientations.
        Sum of each column = 0 (every face has outward normal for one cell, inward for the other).

    EXACTNESS:
        d₂ d₁ = 0  (divergence of curl is zero)
    """
    n_cells = len(cell_face_incidence)
    d2 = np.zeros((n_cells, n_faces))

    for cell_idx, face_list in enumerate(cell_face_incidence):
        for face_idx, orientation in face_list:
            # Use += for robustness against accidental duplicates in cell_face_incidence.
            # The verification below will catch any duplicates via col_nonzero != 2.
            d2[cell_idx, face_idx] += orientation

    # Verify: each face in exactly 2 cells with opposite orientations (sum=0)
    if verify:
        col_nonzero = np.sum(d2 != 0, axis=0)  # count nonzero per column
        col_sums = np.sum(d2, axis=0)          # sum per column

        bad_count = np.where(col_nonzero != 2)[0]
        bad_sum = np.where(np.abs(col_sums) > EPS_CLOSE)[0]

        if len(bad_count) > 0:
            raise ValueError(
                f"d₂ verification failed: {len(bad_count)} faces not in exactly 2 cells. "
                f"First bad face: {bad_count[0]} appears in {int(col_nonzero[bad_count[0]])} cells"
            )
        if len(bad_sum) > 0:
            raise ValueError(
                f"d₂ verification failed: {len(bad_sum)} faces have non-zero orientation sum. "
                f"First bad face: {bad_sum[0]} has sum={col_sums[bad_sum[0]]}"
            )

    return d2


def build_incidence_matrices(vertices: np.ndarray,
                             edges: List[Tuple[int, int]],
                             faces: List[List[int]]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build both incidence matrices d₀ and d₁.

    EXACTNESS THEOREM:
        d₁ d₀ = 0

    This is the discrete analog of "curl(grad) = 0".
    Proof: (d₁d₀)ᵢⱼ sums contributions from edges around face i
           coming from vertex j. Each vertex appears in 0 or 2
           consecutive edges of a face, with opposite signs.

    Args:
        vertices, edges, faces: from build_kelvin_cell()

    Returns:
        d0: (E, V) gradient matrix
        d1: (F, E) curl matrix

    VERIFICATION:
        ||d₁d₀|| should be exactly 0.
    """
    d0 = build_d0(vertices, edges)
    d1 = build_d1(vertices, edges, faces)

    # Verify exactness
    d1d0 = d1 @ d0
    if not np.allclose(d1d0, 0):
        raise ValueError(f"Exactness failed: ||d₁d₀|| = {np.linalg.norm(d1d0)}")

    return d0, d1


def build_hodge_laplacian(d0: np.ndarray,
                          d1: np.ndarray,
                          faces_per_edge: int = 2) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build Hodge Laplacian L₁ on edges.

    DEFINITION:
        L₁ = d₀d₀ᵀ + d₁ᵀd₁

    INTERPRETATION:
        d₀d₀ᵀ: "vertex term" - measures disagreement with neighbors
        d₁ᵀd₁: "face term" - measures curl around faces

    Args:
        d0: (E, V) gradient matrix
        d1: (F, E) curl matrix
        faces_per_edge: from mesh contract (2 for surface, 3 for foam)

    Returns:
        L1: (E, E) Hodge Laplacian
        d0d0t: (E, E) vertex term
        d1td1: (E, E) face term

    TRACE IDENTITIES:
        Tr(d₀d₀ᵀ) = 2E           (ALWAYS - each edge has 2 endpoints)
        Tr(d₁ᵀd₁) = k·E          (k = faces_per_edge)
        Tr(L₁) = (2+k)·E

    PROOF of Tr(d₁ᵀd₁) = k·E:
        (d₁ᵀd₁)ₑₑ = Σ_f d₁[f,e]²
        Each edge appears in exactly k faces,
        with coefficient ±1 in each → contribution = k
        Sum over E edges → k·E
    """
    E = d0.shape[0]
    k = faces_per_edge

    d0d0t = d0 @ d0.T
    d1td1 = d1.T @ d1
    L1 = d0d0t + d1td1

    # Verify trace identities (using if/raise for -O robustness)
    tr_d0d0t = np.trace(d0d0t)
    tr_d1td1 = np.trace(d1td1)
    tr_L1 = np.trace(L1)

    if abs(tr_d0d0t - 2*E) >= EPS_CLOSE:
        raise ValueError(f"Tr(d₀d₀ᵀ) = {tr_d0d0t}, expected 2E = {2*E}")
    if abs(tr_d1td1 - k*E) >= EPS_CLOSE:
        raise ValueError(f"Tr(d₁ᵀd₁) = {tr_d1td1}, expected {k}E = {k*E}")
    if abs(tr_L1 - (2+k)*E) >= EPS_CLOSE:
        raise ValueError(f"Tr(L₁) = {tr_L1}, expected {2+k}E = {(2+k)*E}")

    return L1, d0d0t, d1td1


def count_connected_components(d0: np.ndarray) -> int:
    """
    Count connected components of graph via BFS on d₀.

    Args:
        d0: (E, V) gradient matrix

    Returns:
        c: number of connected components

    NOTE:
        Uses adjacency structure from d₀.
        For Kelvin cell, c = 1 (connected graph).
    """
    E, V = d0.shape

    # Build adjacency list from d0
    adj = [[] for _ in range(V)]
    for e in range(E):
        row = d0[e, :]
        verts = np.where(row != 0)[0]
        if len(verts) == 2:
            i, j = verts
            adj[i].append(j)
            adj[j].append(i)

    # BFS to count components (deque for O(1) popleft)
    visited = [False] * V
    components = 0

    for start in range(V):
        if visited[start]:
            continue
        # BFS from start
        queue = deque([start])
        visited[start] = True
        while queue:
            v = queue.popleft()
            for neighbor in adj[v]:
                if not visited[neighbor]:
                    visited[neighbor] = True
                    queue.append(neighbor)
        components += 1

    return components


def get_cycle_space(d0: np.ndarray, tol: float = None,
                    strict: bool = True) -> np.ndarray:
    """
    Compute cycle space H = ker(d₀ᵀ) of the GRAPH.

    IMPORTANT: This is the GRAPH cycle space, NOT surface homology H₁.
        - H(graph) = ker(d₀ᵀ) = dim 13 for Kelvin
        - H₁(surface) = 0 for Kelvin (sphere topology)

    These are different objects. We compute the graph cycle space.

    DEFINITION:
        H is the space of "closed 1-forms" (divergence-free on vertices).
        Equivalently: edge-weighted cycles in the 1-skeleton graph.

    DIMENSION THEOREM:
        dim(ker(d₀ᵀ)) = E - rank(d₀ᵀ) = E - (V - c) = E - V + c

        where c = number of connected components.
        For Kelvin (connected): c = 1, so dim(H) = 36 - 24 + 1 = 13

    Args:
        d0: (E, V) gradient matrix
        tol: tolerance for zero eigenvalues (defaults to EPS_CLOSE from constants)
        strict: If True (default), verify dimension formula.
                If False, skip verification (for exploration/audit).

    Returns:
        H: (E, dim_H) matrix, columns form orthonormal basis of cycle space

    NOTE:
        We compute ker(d₀ᵀ) via eigendecomposition of d₀d₀ᵀ
        since ker(d₀ᵀ) = ker(d₀d₀ᵀ) for real matrices.
    """
    if tol is None:
        tol = EPS_CLOSE

    E, V = d0.shape

    d0d0t = d0 @ d0.T
    eigvals, eigvecs = np.linalg.eigh(d0d0t)

    # Select zero eigenspace
    H = eigvecs[:, eigvals < tol]

    # Verify dimension using general formula with component count
    if strict:
        c = count_connected_components(d0)
        expected_dim = E - V + c
        actual_dim = H.shape[1]
        if actual_dim != expected_dim:
            raise ValueError(
                f"dim(H) = {actual_dim}, expected E - V + c = {E} - {V} + {c} = {expected_dim}"
            )

    return H


def verify_hodge_on_cycle_space(L1: np.ndarray,
                                d0d0t: np.ndarray,
                                H: np.ndarray) -> Dict[str, float]:
    """
    Verify that d₀d₀ᵀ vanishes on cycle space H.

    THEOREM:
        If h ∈ ker(d₀ᵀ), then d₀d₀ᵀ h = 0.

    This means L₁|_H = (d₁ᵀd₁)|_H  (only the face term survives).

    Args:
        L1: Hodge Laplacian
        d0d0t: vertex term
        H: cycle space basis

    Returns:
        dict with verification results
    """
    d0d0t_on_H = H.T @ d0d0t @ H
    L1_on_H = H.T @ L1 @ H

    norm_d0d0t_on_H = np.linalg.norm(d0d0t_on_H)

    return {
        'd0d0t_on_H_norm': norm_d0d0t_on_H,
        'd0d0t_vanishes': norm_d0d0t_on_H < EPS_CLOSE,
        'L1_on_H_trace': np.trace(L1_on_H)
    }


# =============================================================================
# UNIVERSAL VERIFICATION HELPERS
# =============================================================================

def verify_faces_per_edge(d1: np.ndarray, k: int) -> Dict[str, Any]:
    """
    Universal check: every edge has exactly k incident faces.

    This is THE structural invariant for cell complexes:
        - k=2: 2-manifold surface (each edge bounds 2 faces)
        - k=3: foam/Plateau structure (Plateau border = 3 films meet)
        - k=4: SC solid (4 cubes share each edge)

    Args:
        d1: (F, E) face-edge incidence matrix
        k: expected faces per edge from complex type

    Returns:
        dict with:
            'valid': bool - all edges have exactly k faces
            'min': int - minimum faces per edge
            'max': int - maximum faces per edge
            'expected': int - k
            'histogram': dict - {count: n_edges_with_that_count}

    Raises:
        ValueError: if not all edges have k faces (unless called with raise_on_fail=False)

    USAGE:
        # In any operator code:
        result = verify_faces_per_edge(d1, k)
        if not result['valid']:
            print(f"WARNING: {result['histogram']}")
    """
    faces_per_edge_actual = np.sum(np.abs(d1), axis=0)

    fpe_min = int(np.min(faces_per_edge_actual))
    fpe_max = int(np.max(faces_per_edge_actual))
    is_uniform = (fpe_min == fpe_max == k)

    # Build histogram
    unique, counts = np.unique(faces_per_edge_actual, return_counts=True)
    histogram = {int(u): int(c) for u, c in zip(unique, counts)}

    return {
        'valid': is_uniform,
        'min': fpe_min,
        'max': fpe_max,
        'expected': k,
        'histogram': histogram,
    }


def assert_faces_per_edge(d1: np.ndarray, k: int, context: str = "") -> None:
    """
    Assert that all edges have exactly k incident faces.

    This is the fail-fast version of verify_faces_per_edge.
    Call this in any code that depends on the k-uniformity invariant.

    Args:
        d1: (F, E) face-edge incidence matrix
        k: expected faces per edge
        context: optional context string for error message

    Raises:
        ValueError: if invariant violated
    """
    result = verify_faces_per_edge(d1, k)
    if not result['valid']:
        ctx = f" [{context}]" if context else ""
        raise ValueError(
            f"faces_per_edge invariant violated{ctx}: "
            f"expected all edges to have {k} faces, "
            f"got min={result['min']}, max={result['max']}. "
            f"Histogram: {result['histogram']}"
        )


# =============================================================================
# CONTRACT-AWARE WRAPPER
# =============================================================================

def build_operators_from_mesh(mesh: dict) -> dict:
    """
    Build all DEC operators from a contract-compliant mesh dict.

    Args:
        mesh: Contract-compliant mesh dict with V, E, F, faces_per_edge

    Returns:
        dict with:
            d0, d1: incidence matrices
            L1, d0d0t, d1td1: Hodge Laplacian components
            H: cycle space basis
            traces: dict of trace values
            faces_per_edge_histogram: edge-count verification

    VERIFICATION:
        Checks that all edges have exactly k faces (faces_per_edge from contract).
        This is the structural invariant for the complex type.
    """
    # Use clear names: verts/edges/faces (not V/E/F which are also counts)
    verts = mesh['V']
    edges = mesh['E']
    faces = mesh['F']
    k = mesh['faces_per_edge']

    d0, d1 = build_incidence_matrices(verts, edges, faces)
    L1, d0d0t, d1td1 = build_hodge_laplacian(d0, d1, faces_per_edge=k)
    H = get_cycle_space(d0)

    n_E = len(edges)

    # Verify faces_per_edge using universal helper (fail-fast on violation)
    assert_faces_per_edge(d1, k, context=mesh.get('name', 'mesh'))
    fpe_result = verify_faces_per_edge(d1, k)

    return {
        'd0': d0,
        'd1': d1,
        'L1': L1,
        'd0d0t': d0d0t,
        'd1td1': d1td1,
        'H': H,
        'traces': {
            'Tr_d0d0t': np.trace(d0d0t),
            'Tr_d1td1': np.trace(d1td1),
            'Tr_L1': np.trace(L1),
            'expected_d0d0t': 2 * n_E,
            'expected_d1td1': k * n_E,
            'expected_L1': (2 + k) * n_E,
        },
        'faces_per_edge_histogram': fpe_result,
    }


# Self-test when run directly
# Run with: python -m core_math.operators.incidence (from ST_8/)
def _test_init():
    """Incidence operator tests on Kelvin cell (surface, k=2) + BCC foam (k=3)."""
    from ..builders.kelvin import build_kelvin_cell
    from ..builders.multicell_periodic import build_bcc_supercell_periodic

    print("incidence init tests")
    print("-" * 40)

    # --- Surface (k=2): Kelvin cell ---
    verts, edges, faces, _ = build_kelvin_cell()
    V, E, F = len(verts), len(edges), len(faces)

    d0 = build_d0(verts, edges)
    d1 = build_d1(verts, edges, faces)

    # T1: shapes
    assert d0.shape == (E, V), f"d0 shape {d0.shape} != ({E},{V})"
    assert d1.shape == (F, E), f"d1 shape {d1.shape} != ({F},{E})"
    print(f"  T1 shapes: d0=({E},{V}), d1=({F},{E}) OK")

    # T2: exactness d1 @ d0 = 0
    product = d1 @ d0
    assert np.max(np.abs(product)) < EPS_CLOSE, \
        f"d1@d0 != 0: max={np.max(np.abs(product)):.2e}"
    print(f"  T2 exactness: ||d1@d0||={np.max(np.abs(product)):.2e} OK")

    # T3: trace identities (surface k=2)
    tr_d0d0t = np.trace(d0 @ d0.T)
    tr_d1td1 = np.trace(d1.T @ d1)
    assert abs(tr_d0d0t - 2 * E) < EPS_CLOSE, f"Tr(d0d0T)={tr_d0d0t} != 2E={2*E}"
    assert abs(tr_d1td1 - 2 * E) < EPS_CLOSE, f"Tr(d1Td1)={tr_d1td1} != 2E={2*E}"
    print(f"  T3 traces: Tr(d0d0T)={tr_d0d0t:.0f}=2E, Tr(d1Td1)={tr_d1td1:.0f}=2E OK")

    # T4: d0 entries are {-1, 0, +1}
    d0_vals = set(d0.flatten().astype(int))
    assert d0_vals <= {-1, 0, 1}, f"d0 has values {d0_vals}"
    d1_vals = set(d1.flatten().astype(int))
    assert d1_vals <= {-1, 0, 1}, f"d1 has values {d1_vals}"
    print(f"  T4 entries: d0,d1 ∈ {{-1,0,+1}} OK")

    # T5: cycle space dimension
    dim_H = E - V + 1  # graph cycle space
    n_components = count_connected_components(d0)
    assert n_components == 1, f"Kelvin has {n_components} components"
    print(f"  T5 cycle space: dim(H)={dim_H}, connected OK")

    # T9: Hodge verification — d0d0T vanishes on cycle space
    H_basis = get_cycle_space(d0)
    assert H_basis.shape[1] == dim_H, f"H cols {H_basis.shape[1]} != {dim_H}"
    d0d0t = d0 @ d0.T  # (E, E)
    L1 = build_hodge_laplacian(d0, d1, faces_per_edge=2)
    hodge_result = verify_hodge_on_cycle_space(L1, d0d0t, H_basis)
    assert hodge_result['d0d0t_on_H_norm'] < EPS_CLOSE, \
        f"d0d0T on H: {hodge_result['d0d0t_on_H_norm']:.2e}"
    print(f"  T9 Hodge: ||d0d0T on H||={hodge_result['d0d0t_on_H_norm']:.2e} OK")

    # --- Foam (k=3): BCC N=2 ---
    verts_f, edges_f, faces_f, cfi = build_bcc_supercell_periodic(2)
    V_f, E_f, F_f = len(verts_f), len(edges_f), len(faces_f)

    d0_f = build_d0(verts_f, edges_f)
    d1_f = build_d1(verts_f, edges_f, faces_f)

    # T6: exactness on foam
    product_f = d1_f @ d0_f
    assert np.max(np.abs(product_f)) < EPS_CLOSE, \
        f"Foam d1@d0 != 0: {np.max(np.abs(product_f)):.2e}"
    print(f"  T6 foam exactness: ||d1@d0||={np.max(np.abs(product_f)):.2e} OK")

    # T7: trace identities (foam k=3)
    tr_d0d0t_f = np.trace(d0_f @ d0_f.T)
    tr_d1td1_f = np.trace(d1_f.T @ d1_f)
    assert abs(tr_d0d0t_f - 2 * E_f) < EPS_CLOSE, \
        f"Foam Tr(d0d0T)={tr_d0d0t_f} != 2E={2*E_f}"
    assert abs(tr_d1td1_f - 3 * E_f) < EPS_CLOSE, \
        f"Foam Tr(d1Td1)={tr_d1td1_f} != 3E={3*E_f}"
    print(f"  T7 foam traces: Tr(d0d0T)={tr_d0d0t_f:.0f}=2E, "
          f"Tr(d1Td1)={tr_d1td1_f:.0f}=3E OK")

    # T8: d2 on foam
    d2 = build_d2(cfi, F_f)
    C = len(cfi)
    assert d2.shape == (C, F_f), f"d2 shape {d2.shape} != ({C},{F_f})"
    product_d2d1 = d2 @ d1_f
    assert np.max(np.abs(product_d2d1)) < EPS_CLOSE, \
        f"d2@d1 != 0: {np.max(np.abs(product_d2d1)):.2e}"
    print(f"  T8 d2: shape ({C},{F_f}), ||d2@d1||={np.max(np.abs(product_d2d1)):.2e} OK")

    print("-" * 40)
    print("ALL PASS")


if __name__ == "__main__":
    _test_init()
