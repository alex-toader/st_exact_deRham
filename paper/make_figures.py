"""
Generate Figure 1 for the CR Math note.

Band structure Gamma-X-M-R-Gamma, Kelvin N=2.
(a) Exact DEC: clean bands with twofold acoustic degeneracy.
(b) Standard DEC: spurious bands collapsing toward zero.

Usage:
    python3 make_figures.py

Outputs fig_bands.pdf in paper/ directory.
Requires: matplotlib, numpy, scipy

Source: adapted from st_exact_deRham/paper/make_figures.py (fig3)
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

from physics.hodge import build_kelvin_with_dual_info, build_hodge_stars_voronoi
from physics.gauge_bloch import compute_edge_shifts, build_d0_bloch, build_d1_bloch_exact
from physics.bloch import (build_d1_bloch_standard, compute_edge_crossings,
                            build_edge_lookup)

OUT = os.path.abspath(os.path.dirname(__file__))


def fig_bands():
    """Figure 1: Band structure, exact vs standard (2 panels)."""
    print("Figure 1: Band structure...")

    data = build_kelvin_with_dual_info(N=2)
    V, E, F = data['V'], data['E'], data['F']
    L, L_vec = data['L'], data['L_vec']
    star1, star2 = build_hodge_stars_voronoi(data)
    shifts = compute_edge_shifts(V, E, L_vec)
    crossings = compute_edge_crossings(V, E, L)
    edge_lookup = build_edge_lookup(E, crossings)
    k_scale = 2 * np.pi / L

    # BZ path: Gamma -> X -> M -> R -> Gamma
    path_segments = [
        ('$\\Gamma$', 'X', np.array([0,0,0.]), np.array([0.5,0,0.])),
        ('X', 'M', np.array([0.5,0,0.]), np.array([0.5,0.5,0.])),
        ('M', 'R', np.array([0.5,0.5,0.]), np.array([0.5,0.5,0.5])),
        ('R', '$\\Gamma$', np.array([0.5,0.5,0.5]), np.array([0,0,0.])),
    ]

    n_pts_per_seg = 30
    all_eigs_ex = []
    all_eigs_std = []
    all_t = []
    tick_positions = [0]
    tick_labels = ['$\\Gamma$']
    t_offset = 0

    for _, seg_label_end, frac_start, frac_end in path_segments:
        for i in range(n_pts_per_seg):
            alpha = i / n_pts_per_seg
            frac = frac_start + alpha * (frac_end - frac_start)
            k = k_scale * frac

            # Skip Gamma exactly (k=0 has different kernel)
            if np.linalg.norm(k) < 1e-10:
                k = k_scale * 0.001 * np.array([1, 1, 1])

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

            t = t_offset + alpha * np.linalg.norm(frac_end - frac_start)
            all_t.append(t)

            # First 12 bands above zero
            thresh = max(np.max(np.abs(eigs_ex)) * 1e-12, 1e-10)
            phys_ex = eigs_ex[np.abs(eigs_ex) > thresh][:12]
            phys_std = eigs_std[np.abs(eigs_std) > thresh][:12]
            all_eigs_ex.append(phys_ex)
            all_eigs_std.append(phys_std)

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

    # Normalize by 4pi^2/L^2
    norm_factor = (2 * np.pi / L) ** 2
    eigs_ex_arr /= norm_factor
    eigs_std_arr /= norm_factor

    # Two panels
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True)

    # Panel (a): Exact
    ax = axes[0]
    for b in range(min(8, max_bands)):
        ax.plot(t_arr, eigs_ex_arr[:, b], '-', color='#2166ac', lw=0.8, alpha=0.8)
    ax.set_ylabel('$\\omega^2 \\;/\\; (4\\pi^2/L^2)$', fontsize=11)
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

    fig.tight_layout()
    path = os.path.join(OUT, 'fig_bands.pdf')
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"  Saved: {path}")


if __name__ == '__main__':
    fig_bands()
