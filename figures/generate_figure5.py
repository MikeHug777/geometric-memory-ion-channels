#!/usr/bin/env python3
"""
Figure 5: Three convergent derivation routes
Paper: "Pore Symmetry Determines Fractal Memory in Ion Channels"

Schematic showing three independent paths to H ≈ 0.83 for C4:
Burnside-Bouchaud (combinatorics), Directed Percolation (critical phenomena),
XXZ spin chain (quantum many-body). All converge on [0.80, 0.92].
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif'],
    'font.size': 9,
    'axes.labelsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# ── Colors ─────────────────────────────────────────────────────────────────
COL_BB = '#D73027'    # Burnside-Bouchaud — red
COL_DP = '#2166AC'    # Directed Percolation — blue
COL_XXZ = '#66BD63'   # XXZ spin chain — green
COL_CONV = '#333333'  # Convergence zone

# ── FIGURE ─────────────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 4.0),
                                gridspec_kw={'width_ratios': [1.4, 1]},
                                constrained_layout=True)

# ═══════════════════════════════════════════════════════════════════════════
# Panel (a): Convergence funnel diagram
# ═══════════════════════════════════════════════════════════════════════════
ax1.set_xlim(-0.7, 3.7)
ax1.set_ylim(-0.3, 3.8)
ax1.set_xticks([])
ax1.set_yticks([])
for spine in ax1.spines.values():
    spine.set_visible(False)

# Three source boxes at top
boxes = [
    {'x': 0.0, 'y': 3.2, 'label': 'Burnside-\nBouchaud', 'sub': 'Combinatorics +\nstatistical mechanics',
     'color': COL_BB, 'pred': '$H = 0.833$'},
    {'x': 1.5, 'y': 3.2, 'label': 'Directed\nPercolation', 'sub': 'Non-equilibrium\ncritical phenomena',
     'color': COL_DP, 'pred': '$H = 0.84$–$0.92$'},
    {'x': 3.0, 'y': 3.2, 'label': 'XXZ Spin\nChain', 'sub': 'Quantum\nmany-body theory',
     'color': COL_XXZ, 'pred': '$H = 0.80$–$0.85$'},
]

for box in boxes:
    # Main box
    rect = FancyBboxPatch((box['x'] - 0.55, box['y'] - 0.25), 1.1, 0.55,
                           boxstyle='round,pad=0.05',
                           facecolor=box['color'], alpha=0.15,
                           edgecolor=box['color'], linewidth=1.5)
    ax1.add_patch(rect)
    ax1.text(box['x'], box['y'] + 0.02, box['label'], fontsize=8,
             ha='center', va='center', fontweight='bold', color=box['color'])

    # Subtitle
    ax1.text(box['x'], box['y'] - 0.55, box['sub'], fontsize=6,
             ha='center', va='top', color='#777777', style='italic')

    # Prediction
    ax1.text(box['x'], box['y'] - 1.0, box['pred'], fontsize=8,
             ha='center', va='top', color=box['color'], fontweight='bold')

    # Arrow to convergence zone
    ax1.annotate('', xy=(1.5, 0.80), xytext=(box['x'], box['y'] - 1.15),
                 arrowprops=dict(arrowstyle='->', color=box['color'],
                                 lw=1.5, connectionstyle='arc3,rad=0'),
                 zorder=8)

# Convergence zone
conv_rect = FancyBboxPatch((0.2, 0.0), 2.6, 0.7,
                            boxstyle='round,pad=0.1',
                            facecolor='#FFF3E0', edgecolor='#D73027',
                            linewidth=1.2, zorder=5)
ax1.add_patch(conv_rect)
ax1.text(1.5, 0.45, 'Convergence: $H \\in [0.80, 0.92]$',
         fontsize=10, ha='center', va='center', fontweight='bold',
         color=COL_CONV, zorder=6)
ax1.text(1.5, 0.15, 'for $C_4$ ($B = 6$)',
         fontsize=8, ha='center', va='center', color='#666666', zorder=6)

# Measured data annotation
ax1.text(1.5, -0.25, 'Measured: BK $H = 0.81$–$0.93$ (Wawrzkiewicz-Jalowiecka et al., 2020)',
         fontsize=6.5, ha='center', va='top', color='#888888', style='italic')

ax1.text(-0.02, 1.08, 'a', transform=ax1.transAxes, fontsize=12,
         fontweight='bold', va='top')

# ═══════════════════════════════════════════════════════════════════════════
# Panel (b): H-axis comparison
# ═══════════════════════════════════════════════════════════════════════════

# Vertical H axis showing predictions from each route
ax2.set_xlim(-0.5, 3.5)
ax2.set_ylim(0.70, 1.00)
ax2.set_ylabel('Hurst exponent $H$')
ax2.set_xticks([0, 1, 2])
ax2.set_xticklabels(['Burnside-\nBouchaud', 'Directed\nPercolation', 'XXZ\nSpin Chain'],
                      fontsize=7)

# Measured range (BK data)
ax2.axhspan(0.81, 0.93, alpha=0.12, color='#AAAAAA', zorder=0)
ax2.text(2.8, 0.87, 'BK\ndata', fontsize=6, color='#888888', ha='left',
         va='center', style='italic')

# Burnside-Bouchaud: exact prediction H = 0.833
ax2.plot(0, 0.833, 'o', color=COL_BB, markersize=10,
         markeredgecolor='black', markeredgewidth=1.0, zorder=10)
ax2.plot([0, 0], [0.833, 0.833], '-', color=COL_BB, linewidth=2, zorder=5)

# DP: range 0.84–0.92
ax2.plot(1, 0.843, 'v', color=COL_DP, markersize=8,
         markeredgecolor='black', markeredgewidth=0.8, zorder=10)
ax2.plot(1, 0.920, '^', color=COL_DP, markersize=8,
         markeredgecolor='black', markeredgewidth=0.8, zorder=10)
ax2.plot([1, 1], [0.843, 0.920], '-', color=COL_DP, linewidth=2.5, zorder=5)
ax2.text(1.15, 0.843, '$1 - \\theta/2$', fontsize=6, va='center', color=COL_DP)
ax2.text(1.15, 0.920, '$(3-\\alpha)/2$', fontsize=6, va='center', color=COL_DP)

# XXZ: range 0.80–0.85
ax2.plot(2, 0.80, 'v', color=COL_XXZ, markersize=8,
         markeredgecolor='black', markeredgewidth=0.8, zorder=10)
ax2.plot(2, 0.85, '^', color=COL_XXZ, markersize=8,
         markeredgecolor='black', markeredgewidth=0.8, zorder=10)
ax2.plot([2, 2], [0.80, 0.85], '-', color=COL_XXZ, linewidth=2.5, zorder=5)
ax2.text(2.15, 0.80, '$\\Delta_{\\mathrm{dip}}=-0.47$', fontsize=5.5,
         va='center', color=COL_XXZ)
ax2.text(2.15, 0.85, '$\\Delta_{\\mathrm{dip}}=-0.36$', fontsize=5.5,
         va='center', color=COL_XXZ)

# Burnside prediction reference line
ax2.axhline(0.833, color=COL_BB, linestyle='--', linewidth=0.6, alpha=0.5, zorder=1)
ax2.text(-0.45, 0.845, '$H = 0.833$', fontsize=6, va='bottom', color=COL_BB)

ax2.text(-0.12, 1.08, 'b', transform=ax2.transAxes, fontsize=12,
         fontweight='bold', va='top')

# ── Save ───────────────────────────────────────────────────────────────────
outpath = os.path.join(os.path.dirname(__file__), 'figure5_three_routes.png')
fig.savefig(outpath, dpi=300)
fig.savefig(outpath.replace('.png', '.pdf'))
print(f"Saved: {outpath}")
print(f"Saved: {outpath.replace('.png', '.pdf')}")
plt.close()
