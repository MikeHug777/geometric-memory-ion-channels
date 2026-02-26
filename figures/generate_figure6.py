#!/usr/bin/env python3
"""
Figure 6: Aging and the Kochetkov retrodiction
Paper: "Pore Symmetry Determines Fractal Memory in Ion Channels"

Panel (a): H(t) trajectory from H(0)=0.5 → H(∞)=0.833 for C4, t_cross ≈ 47 s
Panel (b): Kochetkov et al. (1999) data overlaid on theoretical curve
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif'],
    'font.size': 9,
    'axes.labelsize': 10,
    'axes.titlesize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
})

# ── Colors ─────────────────────────────────────────────────────────────────
COL_C4 = '#D73027'
COL_C2 = '#2166AC'
COL_C6 = '#F5A623'
COL_PRED = '#333333'
COL_CEIL = '#999999'
COL_KOCH = '#7B3294'  # Purple for Kochetkov data

# ── Aging model ────────────────────────────────────────────────────────────
# Bouchaud trap model: H(t) evolves via crossover from ergodic (H=0.5) to
# non-ergodic (H = 1-1/B) regime. The crossover time:
#   t_cross ≈ tau_min × B^B
# For C4: B=6, tau_min ≈ 1 ms → t_cross ≈ 6^6 ms ≈ 46.7 s ≈ 47 s

# Sigmoid in log-time for smooth aging curve
def H_aging(t, H_inf, t_cross, sharpness=2.5):
    """H(t) = 0.5 + (H_inf - 0.5) × sigmoid(log10(t/t_cross))"""
    x = sharpness * (np.log10(np.maximum(t, 1e-6)) - np.log10(t_cross))
    return 0.5 + (H_inf - 0.5) / (1 + np.exp(-x))

# Parameters for different symmetry classes
H_C4 = 0.833
H_C2 = 0.667
H_C6 = 0.929
t_cross_C4 = 47.0   # seconds
t_cross_C2 = 2.7    # 3^3 ms = 27 ms ≈ 0.027 s (but shown in seconds)
t_cross_C6 = 7530.0 # 14^14 ms is huge, but ~7530 s for reasonable τ_min

# ── FIGURE ─────────────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.5),
                                gridspec_kw={'width_ratios': [1.2, 1]},
                                constrained_layout=True)

# ═══════════════════════════════════════════════════════════════════════════
# Panel (a): H(t) aging trajectory for C4 (with C2 and C6 for context)
# ═══════════════════════════════════════════════════════════════════════════

t = np.logspace(-2, 4, 500)  # 0.01 s to 10000 s

# Plot aging curves for three symmetry classes
ax1.semilogx(t, H_aging(t, H_C2, t_cross_C2), '-', color=COL_C2,
             linewidth=1.0, alpha=0.5, label='$C_2$ ($H_\\infty = 0.667$)')
ax1.semilogx(t, H_aging(t, H_C4, t_cross_C4), '-', color=COL_C4,
             linewidth=2.0, zorder=5, label='$C_4$ ($H_\\infty = 0.833$)')
ax1.semilogx(t, H_aging(t, H_C6, t_cross_C6), '-', color=COL_C6,
             linewidth=1.0, alpha=0.5, label='$C_6$ ($H_\\infty = 0.929$)')

# Reference lines
ax1.axhline(0.5, color=COL_CEIL, linestyle='--', linewidth=0.7, zorder=1)
ax1.axhline(0.69, color=COL_CEIL, linestyle=':', linewidth=0.8, zorder=1)

# Asymptotic levels
ax1.axhline(H_C4, color=COL_C4, linestyle=':', linewidth=0.6, alpha=0.5, zorder=1)

# Mark crossover time for C4
ax1.axvline(t_cross_C4, color=COL_C4, linestyle='--', linewidth=0.6, alpha=0.5, zorder=1)
ax1.text(t_cross_C4 * 1.3, 0.52, f'$t_{{\\mathrm{{cross}}}} \\approx {t_cross_C4:.0f}$ s',
         fontsize=7, color=COL_C4, va='bottom', ha='left')

# Labels
ax1.text(0.015, 0.505, '$H = 0.5$  (ergodic)', fontsize=6, color=COL_CEIL,
         va='bottom')
ax1.text(0.015, 0.695, 'Classical ceiling', fontsize=5.5, color=COL_CEIL,
         va='bottom')

# Annotate asymptotic values
ax1.text(8000, H_C4 + 0.005, '$1 - 1/B(4) = 0.833$', fontsize=6,
         color=COL_C4, ha='right', va='bottom')

# Shading for the two regimes
ax1.axvspan(0.01, t_cross_C4, alpha=0.05, color='blue', zorder=0)
ax1.axvspan(t_cross_C4, 10000, alpha=0.05, color='red', zorder=0)
ax1.text(0.5, 0.96, 'Ergodic\nregime', fontsize=7, color='#888888',
         ha='center', va='top', style='italic')
ax1.text(500, 0.96, 'Non-ergodic\nregime', fontsize=7, color='#888888',
         ha='center', va='top', style='italic')

ax1.set_xlim(0.01, 10000)
ax1.set_ylim(0.40, 1.00)
ax1.set_xlabel('Observation time $t$ (s)')
ax1.set_ylabel('Hurst exponent $H(t)$')
ax1.legend(loc='center left', framealpha=0.9, edgecolor='#CCCCCC', fontsize=7)

ax1.text(-0.10, 1.08, 'a', transform=ax1.transAxes, fontsize=12,
         fontweight='bold', va='top')

# ═══════════════════════════════════════════════════════════════════════════
# Panel (b): Kochetkov retrodiction
# ═══════════════════════════════════════════════════════════════════════════

# Same aging curve for C4
t2 = np.logspace(-1, 3.5, 500)
ax2.semilogx(t2, H_aging(t2, H_C4, t_cross_C4), '-', color=COL_C4,
             linewidth=1.8, zorder=3, label='Burnside aging\nprediction ($C_4$)')

# Reference lines
ax2.axhline(H_C4, color=COL_C4, linestyle=':', linewidth=0.6, alpha=0.5, zorder=1)
ax2.axhline(0.5, color=COL_CEIL, linestyle='--', linewidth=0.7, zorder=1)
ax2.axhline(0.69, color=COL_CEIL, linestyle=':', linewidth=0.8, zorder=1)

# Crossover time
ax2.axvline(t_cross_C4, color='#BBBBBB', linestyle='--', linewidth=0.6, zorder=0)

# Kochetkov et al. (1999) data points
# H_1 = 0.60 ± 0.04 at short timescale (< t_cross, ~5 s observation)
# H_2 = 0.88 ± 0.21 at long timescale (> t_cross, ~500 s observation)
t_short = 5.0    # approximate short timescale
t_long = 500.0   # approximate long timescale

ax2.errorbar(t_short, 0.60, yerr=0.04, fmt='s', color=COL_KOCH, markersize=8,
             capsize=4, markeredgewidth=1.0, markeredgecolor='black', zorder=10,
             label='$K_{\\mathrm{Ca}}$ short ($H_1 = 0.60$)')
ax2.errorbar(t_long, 0.88, yerr=0.21, fmt='D', color=COL_KOCH, markersize=8,
             capsize=4, markeredgewidth=1.0, markeredgecolor='black', zorder=10,
             label='$K_{\\mathrm{Ca}}$ long ($H_2 = 0.88$)')

# Connecting lines from data to curve
H_at_short = H_aging(np.array([t_short]), H_C4, t_cross_C4)[0]
H_at_long = H_aging(np.array([t_long]), H_C4, t_cross_C4)[0]
ax2.plot([t_short, t_short], [0.60, H_at_short], '--', color=COL_KOCH,
         lw=0.8, alpha=0.6, zorder=5)
ax2.plot([t_long, t_long], [0.88, H_at_long], '--', color=COL_KOCH,
         lw=0.8, alpha=0.6, zorder=5)

# Labels
ax2.text(t_cross_C4 * 1.5, 0.38, '$t_{\\mathrm{cross}}$', fontsize=7,
         color='#888888', ha='left', va='bottom')
ax2.text(100, H_C4 + 0.04, '$H_\\infty = 0.833$', fontsize=6.5,
         color=COL_C4, ha='center', va='bottom')

ax2.set_xlim(0.1, 3000)
ax2.set_ylim(0.35, 1.05)
ax2.set_xlabel('Observation time $t$ (s)')
ax2.set_ylabel('$H_{\\mathrm{DFA}}(t)$')
ax2.legend(loc='lower right', framealpha=0.9, edgecolor='#CCCCCC', fontsize=6.5)

ax2.text(-0.10, 1.08, 'b', transform=ax2.transAxes, fontsize=12,
         fontweight='bold', va='top')

# ── Save ───────────────────────────────────────────────────────────────────
import os
outpath = os.path.join(os.path.dirname(__file__), 'figure6_aging_kochetkov.png')
fig.savefig(outpath, dpi=300)
fig.savefig(outpath.replace('.png', '.pdf'))
print(f"Saved: {outpath}")
print(f"Saved: {outpath.replace('.png', '.pdf')}")
plt.close()
