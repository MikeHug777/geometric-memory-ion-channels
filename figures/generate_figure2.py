#!/usr/bin/env python3
"""
Figure 2: Predictions versus published data
Paper: "Pore Symmetry Determines Fractal Memory in Ion Channels"

Panel (a): H = 1 - 1/B(n) for C2-C6 with published data points
Panel (b): BK voltage series — H vs P_open with DP prediction
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from sympy import totient, divisors

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

# ── Symmetry-class colors (consistent across all figures) ──────────────────
COLORS = {
    2: '#2166AC',   # C2 — blue
    3: '#66BD63',   # C3 — green
    4: '#D73027',   # C4 — red
    5: '#7B3294',   # C5 — purple
    6: '#F5A623',   # C6 — gold
}
GRAY_CLASSICAL = '#BCBCBC'

# ── Burnside formula ───────────────────────────────────────────────────────
def burnside(n):
    """B(n) = (1/n) sum_{d|n} phi(n/d) * 2^d"""
    return sum(int(totient(n // d)) * 2**d for d in divisors(n)) // n

def hurst_pred(n):
    return 1.0 - 1.0 / burnside(n)

# ── Data ───────────────────────────────────────────────────────────────────
ns = [2, 3, 4, 5, 6]
Bs = [burnside(n) for n in ns]
Hs = [hurst_pred(n) for n in ns]

# Published data points
# C2: TREK-2 (WJ-2024)
trek2 = {'n': 2, 'H': 0.66, 'err': 0.05, 'label': 'TREK-2'}

# C4: BK beta4 at three voltages (WJ-2024)
bk_20 = {'n': 4, 'H': 0.93, 'err': 0.03, 'popen': 0.74, 'label': 'BK β4 +20 mV'}
bk_40 = {'n': 4, 'H': 0.81, 'err': 0.07, 'popen': 0.80, 'label': 'BK β4 +40 mV'}
bk_60 = {'n': 4, 'H': 0.80, 'err': 0.07, 'popen': 0.88, 'label': 'BK β4 +60 mV'}

# C4: BK beta2 (WJ-2024), two conditions
bk_b2_ca0  = {'n': 4, 'H': 0.70, 'err': 0.02, 'label': 'BK β2 Ca=0'}
bk_b2_ca10 = {'n': 4, 'H': 0.67, 'err': 0.02, 'label': 'BK β2 Ca=10µM'}

# Own data: C4, clean n=15, mean=0.791, median=0.721, range 0.602-1.019
own_h_values = [0.837, 0.846, 0.748, 0.953, 0.968, 0.718, 0.721, 0.714,
                0.721, 0.687, 0.677, 0.602, 0.641, 1.019, 1.018]

# ── FIGURE ─────────────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.8),
                                gridspec_kw={'width_ratios': [2.2, 1]},
                                constrained_layout=True)

# ═══════════════════════════════════════════════════════════════════════════
# Panel (a): H = 1 - 1/B(n) with data
# ═══════════════════════════════════════════════════════════════════════════

# Shaded zones
ax1.axhspan(0.40, 0.50, color='#E0E0E0', alpha=0.4, zorder=0)
ax1.axhspan(0.50, 0.69, color='#F0E8D0', alpha=0.35, zorder=0)
ax1.text(1.6, 0.45, 'No memory', fontsize=6.5, color='#888888',
         ha='left', va='center', style='italic')
ax1.text(1.6, 0.595, 'Classical\nregime', fontsize=6.5, color='#888888',
         ha='left', va='center', style='italic')

# Classical ceiling line
ax1.axhline(0.69, color='#999999', linestyle=':', linewidth=1.0, zorder=1)
ax1.text(6.75, 0.695, 'Classical ceiling\n(Sigg & Carnevale 2025)',
         fontsize=6, color='#777777', ha='right', va='bottom')

# H = 0.5 line
ax1.axhline(0.50, color='#CCCCCC', linestyle='--', linewidth=0.7, zorder=1)

# Smooth theory curve (interpolated for visual continuity)
n_fine = np.linspace(1.5, 6.5, 200)
# Interpolation through known points using a smooth spline-like curve
from scipy.interpolate import PchipInterpolator
n_pts = np.array([1, 2, 3, 4, 5, 6, 7])
h_pts = np.array([0.5, hurst_pred(2), hurst_pred(3), hurst_pred(4),
                   hurst_pred(5), hurst_pred(6),
                   1 - 1/burnside(7) if False else 0.95])
# Only use real values n=2..6
n_theory = np.array([2, 3, 4, 5, 6])
h_theory = np.array([hurst_pred(n) for n in n_theory])

# Plot theory points connected by line
ax1.plot(n_theory, h_theory, 'k-', linewidth=1.8, zorder=5, label='$H = 1 - 1/B(n)$')

# Theory prediction markers (filled for tested, open for untested)
for n, h in zip(ns, Hs):
    if n in [2, 4]:  # tested
        ax1.plot(n, h, 'ko', markersize=7, zorder=6, markerfacecolor='white',
                 markeredgewidth=1.5)
    else:  # untested
        ax1.plot(n, h, 'ko', markersize=8, zorder=6, markerfacecolor='white',
                 markeredgewidth=1.5, markeredgecolor='#AAAAAA')
        ax1.text(n + 0.15, h + 0.02, '?', fontsize=9, fontweight='bold',
                 color='#999999', ha='left', va='bottom')

# ── Published data points ──────────────────────────────────────────────────

# C2: TREK-2
ax1.errorbar(trek2['n'] - 0.08, trek2['H'], yerr=trek2['err'],
             fmt='o', color=COLORS[2], markersize=6, capsize=3,
             markeredgewidth=0.8, markeredgecolor='black', zorder=10,
             label='TREK-2 (C2)')

# C4: BK beta4 — triangles
for i, bk in enumerate([bk_20, bk_40, bk_60]):
    markers = ['^', 's', 'D']
    lbl = bk['label'] if i == 0 else None
    ax1.errorbar(bk['n'] + 0.08 + i * 0.12, bk['H'], yerr=bk['err'],
                 fmt=markers[i], color=COLORS[4], markersize=5.5, capsize=2.5,
                 markeredgewidth=0.6, markeredgecolor='black', zorder=10,
                 label='BK β4 (C4)' if i == 0 else None)

# C4: BK beta2 — open squares
for i, bk in enumerate([bk_b2_ca0, bk_b2_ca10]):
    ax1.errorbar(bk['n'] - 0.15 - i * 0.12, bk['H'], yerr=bk['err'],
                 fmt='s', color=COLORS[4], markersize=4.5, capsize=2,
                 markeredgewidth=0.6, markeredgecolor='black',
                 markerfacecolor='white', zorder=10,
                 label='BK β2 (C4)' if i == 0 else None)

# Own data as box/whisker at C4
bp = ax1.boxplot([own_h_values], positions=[4.55], widths=0.2,
                  patch_artist=True, zorder=8,
                  boxprops=dict(facecolor=COLORS[4], alpha=0.2, edgecolor=COLORS[4]),
                  whiskerprops=dict(color=COLORS[4], linewidth=0.8),
                  capprops=dict(color=COLORS[4], linewidth=0.8),
                  medianprops=dict(color=COLORS[4], linewidth=1.2),
                  flierprops=dict(marker='.', markersize=3, markerfacecolor=COLORS[4]))
ax1.text(4.55, 0.56, 'Own\ndata', fontsize=6, color=COLORS[4],
         ha='center', va='top', style='italic')

# Axis formatting
ax1.set_xlim(1.5, 6.8)
ax1.set_ylim(0.40, 1.02)
ax1.set_xticks(ns)
ax1.set_xticklabels(['$C_2$\n(dimer)', '$C_3$\n(trimer)', '$C_4$\n(tetramer)',
                      '$C_5$\n(pentamer)', '$C_6$\n(hexamer)'], fontsize=7.5)
ax1.set_ylabel('Hurst exponent  $H$')
ax1.set_xlabel('Pore symmetry')

# Prediction labels on the theory curve
for n, h, b in zip(ns, Hs, Bs):
    if n == 2:
        # Shift right to avoid y-axis
        xt, yt = n + 0.05, h + 0.06
        ha_lbl = 'center'
    elif n == 5:
        # Shift right-down below theory curve to avoid red line and B=14
        xt, yt = n + 0.55, h - 0.06
        ha_lbl = 'center'
    elif n == 6:
        # Shift up-left to avoid overlap with B=8
        xt, yt = n - 0.35, h + 0.045
        ha_lbl = 'center'
    else:
        xt, yt = n - 0.35, h + 0.035
        ha_lbl = 'center'
    ax1.annotate(f'$B={b}$\n$H={h:.3f}$',
                 xy=(n, h), xytext=(xt, yt),
                 fontsize=6, color='#555555', ha=ha_lbl,
                 arrowprops=dict(arrowstyle='-', color='#CCCCCC', lw=0.5))

# Legend
handles = [
    Line2D([0], [0], color='black', linewidth=1.5, label='$H = 1 - 1/B(n)$'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS[2],
           markeredgecolor='black', markersize=6, label='TREK-2 ($C_2$)'),
    Line2D([0], [0], marker='^', color='w', markerfacecolor=COLORS[4],
           markeredgecolor='black', markersize=5.5, label='BK β4 ($C_4$)'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='white',
           markeredgecolor=COLORS[4], markersize=4.5, label='BK β2 ($C_4$)'),
    mpatches.Patch(facecolor=COLORS[4], alpha=0.2, edgecolor=COLORS[4],
                   label='Own data ($C_4$, n=15)'),
]
ax1.legend(handles=handles, loc='lower right', framealpha=0.9, edgecolor='#CCCCCC')

ax1.text(-0.08, 1.06, 'a', transform=ax1.transAxes, fontsize=12,
         fontweight='bold', va='top')

# ═══════════════════════════════════════════════════════════════════════════
# Panel (b): BK voltage series — H vs P_open
# ═══════════════════════════════════════════════════════════════════════════

# DP critical prediction curve
popen_range = np.linspace(0.50, 0.98, 100)
# Simple model: H peaks near P_c ~ 0.7, symmetric Gaussian-like envelope
P_c = 0.70
sigma = 0.15
H_dp_envelope = 0.833 + (0.92 - 0.833) * np.exp(-0.5 * ((popen_range - P_c) / sigma)**2)

ax2.fill_between(popen_range, 0.833, H_dp_envelope, alpha=0.12, color=COLORS[4],
                 zorder=1)
ax2.plot(popen_range, H_dp_envelope, '--', color=COLORS[4], linewidth=1.0,
         alpha=0.6, zorder=2, label='DP critical envelope')

# Burnside asymptote
ax2.axhline(0.833, color='black', linestyle='-.', linewidth=0.8, zorder=3)
ax2.text(0.56, 0.838, '$H_{\\mathrm{Burnside}} = 0.833$', fontsize=6.5,
         ha='left', va='bottom', color='#444444')

# Classical ceiling
ax2.axhline(0.69, color='#999999', linestyle=':', linewidth=0.8, zorder=1)

# Data points: BK beta4
for bk, mk, voltage in [(bk_20, '^', '+20 mV'), (bk_40, 's', '+40 mV'),
                          (bk_60, 'D', '+60 mV')]:
    ax2.errorbar(bk['popen'], bk['H'], yerr=bk['err'],
                 fmt=mk, color=COLORS[4], markersize=6, capsize=3,
                 markeredgewidth=0.7, markeredgecolor='black', zorder=10)
    # Position voltage labels to avoid overlap with Burnside line and each other
    if voltage == '+20 mV':
        tx, ty = bk['popen'] + 0.02, bk['H'] + 0.02
        ha_align = 'left'
    elif voltage == '+40 mV':
        tx, ty = bk['popen'] - 0.06, bk['H'] - 0.025
        ha_align = 'right'
    else:  # +60 mV
        tx, ty = bk['popen'] + 0.02, bk['H'] - 0.025
        ha_align = 'left'
    ax2.annotate(voltage, xy=(bk['popen'], bk['H']),
                 xytext=(tx, ty), fontsize=6.5, ha=ha_align, color='#444444')

# Mark P_c
ax2.axvline(P_c, color='#AAAAAA', linestyle='--', linewidth=0.6, zorder=0)
ax2.text(P_c, 0.63, '$P_c$', fontsize=7, ha='center', color='#888888')

ax2.set_xlim(0.55, 0.98)
ax2.set_ylim(0.62, 1.00)
ax2.set_xlabel('Open probability  $P_{\\mathrm{open}}$')
ax2.set_ylabel('$H_{\\mathrm{DFA}}$')

ax2.text(-0.15, 1.06, 'b', transform=ax2.transAxes, fontsize=12,
         fontweight='bold', va='top')

# ── Save ───────────────────────────────────────────────────────────────────
outpath = '/Users/mikehug/Documents/personality-overview-core/pcm_model/data/Forschung/Interface_First/paper/figures/figure2_predictions_vs_data.png'
fig.savefig(outpath, dpi=300)
outpath_pdf = outpath.replace('.png', '.pdf')
fig.savefig(outpath_pdf)
print(f"Saved: {outpath}")
print(f"Saved: {outpath_pdf}")
plt.close()
