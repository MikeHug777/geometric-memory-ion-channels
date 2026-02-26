#!/usr/bin/env python3
"""
Figure 7: Information-theoretic profiles of pore symmetry classes
Paper: "Pore Symmetry Determines Fractal Memory in Ion Channels"

Bandwidth h (bits/event) vs memory depth D (events) for C2–C6,
showing h × log D → 2. Functional roles annotated.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

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
# Entropy rate: h = 2/(B+2) bits/event
h_vals = [2.0 / (B + 2) for B in Bs]
# Memory depth: D = 10^B events
D_vals = [10**B for B in Bs]
# Bandwidth-storage product: h × log10(D) = 2B/(B+2)
products = [h * np.log10(D) for h, D in zip(h_vals, D_vals)]

# Channel examples
examples = {
    2: 'Hv1, TREK-2',
    3: 'ASIC1a, P2X7',
    4: 'BK, Kv, Nav',
    5: r'$\alpha$7 nAChR',
    6: 'Cx36, Orai',
}

# Functional roles
roles = {
    2: 'Sensor',
    3: 'Sensor',
    4: 'Computer',
    5: 'Transmitter',
    6: 'Transmitter',
}

# ── FIGURE ─────────────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.5),
                                gridspec_kw={'width_ratios': [1.3, 1]},
                                constrained_layout=True)

# ═══════════════════════════════════════════════════════════════════════════
# Panel (a): h vs D (log scale) — the information landscape
# ═══════════════════════════════════════════════════════════════════════════

# Constant product curves (h × log10(D) = const)
D_range = np.logspace(2, 16, 200)
for level in [0.5, 1.0, 1.5, 2.0]:
    h_curve = level / np.log10(D_range)
    valid = (h_curve > 0) & (h_curve < 0.6)
    ax1.plot(D_range[valid], h_curve[valid], '-', color='#DDDDDD',
             linewidth=0.5, zorder=0)

# h × log10(D) = 2 limit curve (prominent)
h_limit = 2.0 / np.log10(D_range)
valid = (h_limit > 0) & (h_limit < 0.6)
ax1.plot(D_range[valid], h_limit[valid], '--', color='#AAAAAA',
         linewidth=1.0, zorder=1, label='$h \\cdot \\log_{10} D = 2$')

# Data points
for n, B, h, D in zip(ns, Bs, h_vals, D_vals):
    ax1.plot(D, h, 'o', color=COLORS[n], markersize=10,
             markeredgewidth=1.2, markeredgecolor='black', zorder=10)

# Annotations with channel examples — clean labels without connector lines
label_pos = {
    2: {'dx': 0.6, 'dy': 0.0, 'ha': 'left'},
    3: {'dx': 0.6, 'dy': 0.0, 'ha': 'left'},
    4: {'dx': 0.6, 'dy': 0.0, 'ha': 'left'},
    5: {'dx': 0.6, 'dy': 0.0, 'ha': 'left'},
    6: {'dx': -0.8, 'dy': 0.0, 'ha': 'right'},
}
for n, B, h, D in zip(ns, Bs, h_vals, D_vals):
    pos = label_pos[n]
    ax1.text(D * 10**pos['dx'], h + pos['dy'],
             f'$C_{n}$  {examples[n]}',
             fontsize=6, color=COLORS[n], ha=pos['ha'], va='center')

# Role regions (shaded)
# Sensor: D < 10^5
ax1.axvspan(1e2, 1e5, alpha=0.06, color='#2166AC', zorder=0)
ax1.text(10**3.5, 0.48, 'Sensor', fontsize=8, color='#2166AC',
         ha='center', va='center', style='italic', fontweight='bold', alpha=0.6)

# Computer: D ~ 10^6
ax1.axvspan(1e5, 1e8, alpha=0.06, color='#D73027', zorder=0)
ax1.text(10**6.5, 0.48, 'Computer', fontsize=8, color='#D73027',
         ha='center', va='center', style='italic', fontweight='bold', alpha=0.6)

# Transmitter: D > 10^8
ax1.axvspan(1e8, 1e16, alpha=0.06, color='#7B3294', zorder=0)
ax1.text(1e11, 0.48, 'Transmitter', fontsize=8, color='#7B3294',
         ha='center', va='center', style='italic', fontweight='bold', alpha=0.6)

ax1.set_xscale('log')
ax1.set_xlim(1e2, 1e16)
ax1.set_ylim(0.05, 0.50)
ax1.set_xlabel('Memory depth $D$ (events)')
ax1.set_ylabel('Entropy rate $h$ (bits/event)')

ax1.text(-0.10, 1.08, 'a', transform=ax1.transAxes, fontsize=12,
         fontweight='bold', va='top')

# ═══════════════════════════════════════════════════════════════════════════
# Panel (b): Bandwidth-storage product h × log10(D) vs n
# ═══════════════════════════════════════════════════════════════════════════

ax2.axhline(2.0, color='#AAAAAA', linestyle='--', linewidth=0.8, zorder=1)
ax2.text(6.0, 2.05, 'Limit = 2', fontsize=6.5, color='#888888',
         ha='right', va='bottom')

for n, B, prod in zip(ns, Bs, products):
    ax2.bar(n, prod, width=0.6, color=COLORS[n], edgecolor='black',
            linewidth=0.6, zorder=5, alpha=0.85)
    ax2.text(n, prod + 0.03, f'{prod:.2f}', fontsize=7, ha='center',
             va='bottom', color='#444444')

ax2.set_xticks(ns)
ax2.set_xticklabels(['$C_2$', '$C_3$', '$C_4$', '$C_5$', '$C_6$'])
ax2.set_xlabel('Pore symmetry')
ax2.set_ylabel('$h \\times \\log_{10} D$')
ax2.set_ylim(0, 2.3)
ax2.set_xlim(1.3, 6.7)

ax2.text(-0.12, 1.08, 'b', transform=ax2.transAxes, fontsize=12,
         fontweight='bold', va='top')

# ── Save ───────────────────────────────────────────────────────────────────
outpath = os.path.join(os.path.dirname(__file__), 'figure7_information_profiles.png')
fig.savefig(outpath, dpi=300)
fig.savefig(outpath.replace('.png', '.pdf'))
print(f"Saved: {outpath}")
print(f"Saved: {outpath.replace('.png', '.pdf')}")
plt.close()
