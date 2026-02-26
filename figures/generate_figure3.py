#!/usr/bin/env python3
"""
Figure 3: Six observable signatures from a single B(n)
Paper: "Pore Symmetry Determines Fractal Memory in Ion Channels"

3×2 grid: H, β, α_F, α_A, R_aging, Δα vs C2–C6
All derive from μ = 2/B(n)
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif'],
    'font.size': 9,
    'axes.labelsize': 9,
    'axes.titlesize': 9,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
})

# ── Colors ─────────────────────────────────────────────────────────────────
COLORS = {
    2: '#2166AC',   # C2 — blue
    3: '#66BD63',   # C3 — green
    4: '#D73027',   # C4 — red
    5: '#7B3294',   # C5 — purple
    6: '#F5A623',   # C6 — gold
}

# ── Data ───────────────────────────────────────────────────────────────────
ns = [2, 3, 4, 5, 6]
Bs = [3, 4, 6, 8, 14]

# Six observables from μ = 2/B(n)
H_vals     = [1 - 1.0/B for B in Bs]           # H = 1 - 1/B
beta_vals  = [1 - 2.0/B for B in Bs]           # β = 1 - 2/B
alphaF_vals = [1 - 2.0/B for B in Bs]          # α_F = 1 - 2/B (same as β)
alphaA_vals = [-2.0/B for B in Bs]             # α_A = -2/B
R_aging_vals = [6, 10, 15, 22, 33]  # Table 2 values: R_aging for 1 h recording
dalpha_vals = [B/2.0 - 1 for B in Bs]          # Δα = B/2 - 1

# ── FIGURE ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), constrained_layout=True)

observables = [
    (H_vals, '$H$', 'Hurst exponent', (0, 1.0), 'a'),
    (beta_vals, r'$\beta$', 'Spectral slope', (0, 1.0), 'b'),
    (alphaF_vals, r'$\alpha_F$', 'Fano exponent', (0, 1.0), 'c'),
    (alphaA_vals, r'$\alpha_A$', 'Allan variance slope', (-0.8, 0), 'd'),
    (R_aging_vals, '$R_{\\mathrm{aging}}$', 'Aging ratio (1 h)', (0, 40), 'e'),
    (dalpha_vals, r'$\Delta\alpha$', 'Multifractal width', (0, 7), 'f'),
]

for idx, (vals, ylabel, title, ylims, label) in enumerate(observables):
    ax = axes[idx // 3, idx % 3]

    # Bar plot
    for i, (n, v) in enumerate(zip(ns, vals)):
        ax.bar(n, v, width=0.6, color=COLORS[n], edgecolor='black',
               linewidth=0.5, alpha=0.85, zorder=5)

    # Value labels on bars
    for n, v in zip(ns, vals):
        offset = 0.02 * (ylims[1] - ylims[0])
        va_pos = 'bottom'
        if v < 0:
            offset = -offset
            va_pos = 'top'
        if title == 'Aging ratio (1 h)':
            ax.text(n, v + offset, f'{v:.0f}x', fontsize=5.5,
                    ha='center', va=va_pos, color='#444444')
        else:
            ax.text(n, v + offset, f'{v:.3f}', fontsize=5.5,
                    ha='center', va=va_pos, color='#444444')

    ax.set_xticks(ns)
    ax.set_xticklabels(['$C_2$', '$C_3$', '$C_4$', '$C_5$', '$C_6$'])
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=8, style='italic', pad=3)
    ax.set_ylim(ylims)
    ax.set_xlim(1.3, 6.7)

    # Reference lines where appropriate
    if ylabel == '$H$':
        ax.axhline(0.5, color='#CCCCCC', linestyle='--', linewidth=0.6)
        ax.axhline(0.69, color='#999999', linestyle=':', linewidth=0.6)
    elif ylabel in [r'$\beta$', r'$\alpha_F$']:
        ax.axhline(0, color='#CCCCCC', linestyle='-', linewidth=0.5)

    # Panel label
    ax.text(-0.12, 1.14, label, transform=ax.transAxes, fontsize=11,
            fontweight='bold', va='top')

# Common formula annotation at bottom
fig.text(0.5, -0.02, 'All six observables derive from $\\mu = 2/B(n)$',
         fontsize=9, ha='center', va='top', style='italic', color='#666666')

# ── Save ───────────────────────────────────────────────────────────────────
outpath = os.path.join(os.path.dirname(__file__), 'figure3_six_observables.png')
fig.savefig(outpath, dpi=300)
fig.savefig(outpath.replace('.png', '.pdf'))
print(f"Saved: {outpath}")
print(f"Saved: {outpath.replace('.png', '.pdf')}")
plt.close()
