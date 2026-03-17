"""
Generate figures for the exact de Rham paper.

Usage:
    python3 make_figures.py           # all figures
    python3 make_figures.py fig2      # only figure 2

Outputs PDF files in paper/ directory.
Requires: matplotlib, numpy, scipy

Mar 2026
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
from scipy.linalg import eigh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from physics.hodge import (build_kelvin_with_dual_info, build_hodge_stars_voronoi)
from physics.gauge_bloch import compute_edge_shifts, build_d0_bloch, build_d1_bloch_exact
from physics.bloch import (build_d1_bloch_standard, compute_edge_crossings,
                            build_edge_lookup)

OUT = os.path.abspath(os.path.dirname(__file__))


def _solve(data, frac, direction, use_exact=True):
    """Solve curl-curl eigenvalue problem at given k."""
    V, E, F = data['V'], data['E'], data['F']
    L, L_vec = data['L'], data['L_vec']
    star1, star2 = build_hodge_stars_voronoi(data)
    shifts = compute_edge_shifts(V, E, L_vec)

    d_norm = np.array(direction, dtype=float)
    d_norm /= np.linalg.norm(d_norm)
    k = 2 * np.pi / L * frac * d_norm
    k2 = np.dot(k, k)

    d0_k = build_d0_bloch(V, E, k, L_vec, shifts)

    if use_exact:
        d1_k = build_d1_bloch_exact(V, E, F, k, L_vec, d0_k)
    else:
        crossings = compute_edge_crossings(V, E, L)
        edge_lookup = build_edge_lookup(E, crossings)
        d1_k = build_d1_bloch_standard(V, E, F, L, k, edge_lookup, crossings)

    K = d1_k.conj().T @ np.diag(star2) @ d1_k
    K = 0.5 * (K + K.conj().T)
    M = np.diag(star1)
    eigs = np.sort(np.real(eigh(K, M, eigvals_only=True)))
    norm_d1d0 = np.linalg.norm(d1_k @ d0_k)

    return eigs, k2, norm_d1d0, len(V)


# ====================================================================
def fig2_convergence():
    """Fig 2: Convergence of c² with N — exact vs standard.

    Log-log plot of |c²-1| vs N for Kelvin (exact: p=2.00, standard: oscillating).
    """
    print("Fig 2: Convergence...")

    frac = 0.05
    Ns_ex = [2, 3, 4, 5]
    Ns_std = [2, 3, 4, 5]
    errors_ex = []
    c2_std_vals = []

    for N in Ns_std:
        data = build_kelvin_with_dual_info(N=N)
        star1, star2 = build_hodge_stars_voronoi(data)

        # Exact (only for Ns_ex)
        if N in Ns_ex:
            eigs_ex, k2, _, nV = _solve(data, frac, [1, 0, 0], use_exact=True)
            thresh = max(np.max(np.abs(eigs_ex)) * 1e-12, 1e-10)
            phys_ex = eigs_ex[np.abs(eigs_ex) > thresh]
            c2_ex = phys_ex[0] / k2
            errors_ex.append(abs(c2_ex - 1.0))
        else:
            eigs_ex, k2, _, nV = _solve(data, frac, [1, 0, 0], use_exact=True)

        # Standard
        eigs_std, _, _, _ = _solve(data, frac, [1, 0, 0], use_exact=False)
        thresh_s = max(np.max(np.abs(eigs_std)) * 1e-12, 1e-10)
        phys_std = eigs_std[np.abs(eigs_std) > thresh_s]
        c2_std = phys_std[0] / k2
        c2_std_vals.append(abs(c2_std - 1.0))

    fig, ax = plt.subplots(1, 1, figsize=(5, 3.5))

    ax.loglog(Ns_ex, errors_ex, 'o-', color='#2166ac', lw=2, ms=7,
              label='Exact DEC')
    ax.loglog(Ns_std, c2_std_vals, 's--', color='#b2182b', lw=2, ms=7,
              label='Standard DEC')

    # Reference slope p=2
    N_ref = np.array([2, 5])
    slope_ref = errors_ex[0] * (N_ref / Ns_ex[0]) ** (-2)
    ax.loglog(N_ref, slope_ref, ':', color='gray', lw=1, label='$O(N^{-2})$')

    # Annotate N=5 standard as spurious
    ax.annotate('spurious\nmode', xy=(5, c2_std_vals[3]),
                xytext=(4.3, 5e-3), fontsize=7, color='#b2182b',
                arrowprops=dict(arrowstyle='->', color='#b2182b', lw=0.8))

    ax.set_xlabel('Supercell size $N$', fontsize=11)
    ax.set_ylabel('$|c^2 - 1|$', fontsize=11)
    ax.set_title('Dispersion convergence (Kelvin, 5% BZ)', fontsize=11)
    ax.legend(fontsize=9, loc='upper right')
    ax.set_xticks(Ns_std)
    ax.set_xticklabels([str(n) for n in Ns_std])
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    path = os.path.join(OUT, 'fig2_convergence.pdf')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {path}")


# ====================================================================
def fig3_band_structure():
    """Fig 3: Band structure Γ→X→M→R→Γ, exact vs standard + ||d1d0||.

    Three panels:
    (a) Exact: clean bands
    (b) Standard: spurious bands
    (c) ||d1d0|| along BZ path
    """
    print("Fig 3: Band structure...")

    data = build_kelvin_with_dual_info(N=2)
    V, E, F = data['V'], data['E'], data['F']
    L, L_vec = data['L'], data['L_vec']
    nV, nE = len(V), len(E)
    star1, star2 = build_hodge_stars_voronoi(data)
    shifts = compute_edge_shifts(V, E, L_vec)
    crossings = compute_edge_crossings(V, E, L)
    edge_lookup = build_edge_lookup(E, crossings)
    k_scale = 2 * np.pi / L

    # BZ path: Γ→X→M→R→Γ
    path_segments = [
        ('$\\Gamma$', 'X', np.array([0,0,0.]), np.array([0.5,0,0.])),
        ('X', 'M', np.array([0.5,0,0.]), np.array([0.5,0.5,0.])),
        ('M', 'R', np.array([0.5,0.5,0.]), np.array([0.5,0.5,0.5])),
        ('R', '$\\Gamma$', np.array([0.5,0.5,0.5]), np.array([0,0,0.])),
    ]

    n_pts_per_seg = 30
    all_eigs_ex = []
    all_eigs_std = []
    all_norms = []
    all_t = []
    tick_positions = [0]
    tick_labels = ['$\\Gamma$']
    t_offset = 0

    for seg_label_start, seg_label_end, frac_start, frac_end in path_segments:
        for i in range(n_pts_per_seg):
            alpha = i / n_pts_per_seg
            frac = frac_start + alpha * (frac_end - frac_start)
            k = k_scale * frac

            # Skip Γ exactly (k=0 has different kernel)
            if np.linalg.norm(k) < 1e-10:
                k = k_scale * 0.001 * np.array([1, 1, 1])

            k2 = np.dot(k, k)
            d0_k = build_d0_bloch(V, E, k, L_vec, shifts)

            # Exact
            d1_ex = build_d1_bloch_exact(V, E, F, k, L_vec, d0_k)
            K_ex = d1_ex.conj().T @ np.diag(star2) @ d1_ex
            K_ex = 0.5 * (K_ex + K_ex.conj().T)
            M = np.diag(star1)
            eigs_ex = np.sort(np.real(eigh(K_ex, M, eigvals_only=True)))

            # Standard
            d1_std = build_d1_bloch_standard(V, E, F, L, k, edge_lookup, crossings)
            K_std = d1_std.conj().T @ np.diag(star2) @ d1_std
            K_std = 0.5 * (K_std + K_std.conj().T)
            eigs_std = np.sort(np.real(eigh(K_std, M, eigvals_only=True)))

            norm_d1d0 = np.linalg.norm(d1_std @ d0_k)

            t = t_offset + alpha * np.linalg.norm(frac_end - frac_start)
            all_t.append(t)
            # Take first 12 bands above zero
            thresh = max(np.max(np.abs(eigs_ex)) * 1e-12, 1e-10)
            phys_ex = eigs_ex[np.abs(eigs_ex) > thresh][:12]
            phys_std = eigs_std[np.abs(eigs_std) > thresh][:12]
            all_eigs_ex.append(phys_ex)
            all_eigs_std.append(phys_std)
            all_norms.append(norm_d1d0)

        t_offset += np.linalg.norm(frac_end - frac_start)
        tick_positions.append(t_offset)
        tick_labels.append(seg_label_end)

    # Pad arrays to equal length
    max_bands = max(len(e) for e in all_eigs_ex)
    def pad(arr_list, max_len):
        result = np.full((len(arr_list), max_len), np.nan)
        for i, arr in enumerate(arr_list):
            result[i, :len(arr)] = arr
        return result

    eigs_ex_arr = pad(all_eigs_ex, max_bands)
    eigs_std_arr = pad(all_eigs_std, max_bands)
    t_arr = np.array(all_t)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=False)

    # Panel (a): Exact
    ax = axes[0]
    for b in range(min(8, max_bands)):
        ax.plot(t_arr, eigs_ex_arr[:, b], '-', color='#2166ac', lw=0.8, alpha=0.8)
    ax.set_ylabel('$\\omega^2$', fontsize=11)
    ax.set_title('(a) Exact DEC', fontsize=11)
    for tp in tick_positions:
        ax.axvline(tp, color='gray', lw=0.3, alpha=0.5)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_xlim(t_arr[0], t_arr[-1])
    ax.set_ylim(bottom=0)

    # Panel (b): Standard
    ax = axes[1]
    for b in range(min(8, max_bands)):
        ax.plot(t_arr, eigs_std_arr[:, b], '-', color='#b2182b', lw=0.8, alpha=0.8)
    ax.set_title('(b) Standard DEC', fontsize=11)
    for tp in tick_positions:
        ax.axvline(tp, color='gray', lw=0.3, alpha=0.5)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_xlim(t_arr[0], t_arr[-1])
    ax.set_ylim(bottom=0)

    # Panel (c): ||d1d0||
    ax = axes[2]
    ax.plot(t_arr, all_norms, '-', color='#762a83', lw=1.5)
    ax.set_ylabel('$\\|d_1^{\\mathrm{std}} d_0\\|$', fontsize=11)
    ax.set_title('(c) Exactness violation', fontsize=11)
    for tp in tick_positions:
        ax.axvline(tp, color='gray', lw=0.3, alpha=0.5)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_xlim(t_arr[0], t_arr[-1])

    fig.tight_layout()
    path = os.path.join(OUT, 'fig3_bands.pdf')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {path}")


# ====================================================================
if __name__ == '__main__':
    targets = sys.argv[1:] if len(sys.argv) > 1 else ['fig2']

    dispatch = {
        'fig2': fig2_convergence,
        'fig3': fig3_band_structure,
    }

    for t in targets:
        if t in dispatch:
            dispatch[t]()
        else:
            print(f"Unknown figure: {t}. Available: {list(dispatch.keys())}")
