#!/usr/bin/env python3
"""
Figure 4: Own patch-clamp data — DFA results
Paper: "Pore Symmetry Determines Fractal Memory in Ion Channels"

Panel (a): Representative DFA log-log plot (HEK293_Cell01 @ 60 mV, H=0.837)
Panel (b): Distribution of H values (n=15) — strip plot with Burnside prediction
Panel (c): H versus recording duration — convergence toward prediction
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys, os

# Add calculations directory to path for DFA function
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'calculations'))
from dfa_hurst_analysis import dfa

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

# ── Colors (consistent with Figure 2) ─────────────────────────────────────
COL_C4 = '#D73027'       # C4 red
COL_PRED = '#333333'     # Burnside prediction
COL_CEIL = '#999999'     # Classical ceiling
COL_HEK = '#D73027'      # HEK293 — red
COL_MOUSE = '#2166AC'    # Mouse CM — blue
COL_RABBIT = '#66BD63'   # Rabbit CM — green

# ── Data directory ─────────────────────────────────────────────────────────
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..',
                        'calculations', 'real_results_det_besttrans')

# ── Parse summary CSV ──────────────────────────────────────────────────────
import csv
records = []
with open(os.path.join(DATA_DIR, 'summary_with_qc.csv')) as f:
    reader = csv.DictReader(f)
    for row in reader:
        name = row['channel_name']
        h = float(row['calculated_h'])
        dev = float(row['deviation_percent'])
        trans = int(row['transitions'])
        popen = float(row['open_fraction'])
        idealized_file = row['idealized_file']

        # Parse cell type
        if name.startswith('HEK'):
            cell_type = 'HEK293'
        elif name.startswith('Mouse'):
            cell_type = 'Mouse CM'
        else:
            cell_type = 'Rabbit CM'

        # Parse duration from name (e.g., "HEK293_Cell01_60mV_160s_raw_DET")
        parts = name.split('_')
        duration = None
        for p in parts:
            if p.endswith('s') and p[:-1].isdigit():
                duration = int(p[:-1])
                break

        # Parse voltage
        voltage = None
        for p in parts:
            if 'mV' in p:
                voltage = int(p.replace('mV', ''))
                break

        records.append({
            'name': name,
            'h': h,
            'dev': dev,
            'trans': trans,
            'popen': popen,
            'cell_type': cell_type,
            'duration': duration,
            'voltage': voltage,
            'idealized_file': idealized_file,
        })

# Exclude outlier (MouseCM_Cell10, H=1.271, transitions=21)
clean = [r for r in records if r['h'] < 1.1]
outlier = [r for r in records if r['h'] >= 1.1]

# ── FIGURE ─────────────────────────────────────────────────────────────────
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(7.2, 3.4),
                                     gridspec_kw={'width_ratios': [1, 0.85, 1.1]},
                                     constrained_layout=True)

# ═══════════════════════════════════════════════════════════════════════════
# Panel (a): DFA log-log plot for best recording
# ═══════════════════════════════════════════════════════════════════════════

# Load best recording: HEK293_Cell01_60mV (H=0.837, dev=0.4%)
best_file = os.path.join(DATA_DIR, 'HEK293_Cell01_60mV_160s_raw_DET_col72_idealized.csv')
best_data = np.loadtxt(best_file)

# Run DFA
scales, fluctuations, H_calc = dfa(best_data)

# Plot DFA points
ax1.loglog(scales, fluctuations, 'o', color=COL_C4, markersize=4,
           markeredgewidth=0.5, markeredgecolor='black', alpha=0.8, zorder=5)

# Fit line
log_s = np.log10(scales)
log_f = np.log10(fluctuations)
p = np.polyfit(log_s, log_f, 1)
fit_x = np.logspace(log_s[0], log_s[-1], 100)
fit_y = 10**(np.polyval(p, np.log10(fit_x)))
ax1.loglog(fit_x, fit_y, '-', color='black', linewidth=1.5, zorder=4,
           label=f'$H = {H_calc:.3f}$')

# Reference slopes
ref_x = np.array([scales[0], scales[-1]])
# H=0.5 reference
y05 = fluctuations[0] * (ref_x / scales[0])**0.5
ax1.loglog(ref_x, y05 * 0.6, '--', color=COL_CEIL, linewidth=0.7, alpha=0.7)
# Place H=0.5 label along the reference line (not at right edge)
mid_x = np.sqrt(ref_x[0] * ref_x[-1])  # geometric midpoint in log space
mid_y05 = fluctuations[0] * (mid_x / scales[0])**0.5 * 0.6
ax1.text(mid_x, mid_y05 * 1.15, '$H=0.5$', fontsize=6,
         color=COL_CEIL, va='bottom', ha='center')

ax1.set_xlabel('Scale (samples)')
ax1.set_ylabel('Fluctuation $F(s)$')
ax1.legend(loc='upper left', fontsize=7, framealpha=0.9)
ax1.set_title('HEK293 Cell01, 60 mV', fontsize=8, style='italic')
ax1.text(-0.12, 1.14, 'a', transform=ax1.transAxes, fontsize=12,
         fontweight='bold', va='top')

# ═══════════════════════════════════════════════════════════════════════════
# Panel (b): Strip plot of H values
# ═══════════════════════════════════════════════════════════════════════════

# Burnside prediction and classical ceiling
ax2.axhline(0.833, color=COL_PRED, linestyle='--', linewidth=1.0, zorder=1)
ax2.axhline(0.69, color=COL_CEIL, linestyle=':', linewidth=0.8, zorder=1)

# Text labels for reference lines
ax2.text(0.98, 0.84, '$H_{\\mathrm{pred}}$', fontsize=6, color=COL_PRED,
         va='bottom', ha='right', transform=ax2.get_yaxis_transform())
ax2.text(0.98, 0.695, 'Class. ceiling', fontsize=5, color=COL_CEIL,
         va='bottom', ha='right', transform=ax2.get_yaxis_transform())

# Jitter the x-position for visibility
np.random.seed(42)
for r in clean:
    if r['cell_type'] == 'HEK293':
        color = COL_HEK
        xbase = 0.0
    elif r['cell_type'] == 'Mouse CM':
        color = COL_MOUSE
        xbase = 0.0
    else:
        color = COL_RABBIT
        xbase = 0.0

    jitter = np.random.uniform(-0.15, 0.15)
    marker = 'o'
    msize = 5

    # Highlight best matches
    if r['dev'] < 2.0:
        edgewidth = 1.5
        msize = 7
    else:
        edgewidth = 0.5

    ax2.plot(xbase + jitter, r['h'], marker, color=color, markersize=msize,
             markeredgewidth=edgewidth, markeredgecolor='black', zorder=10, alpha=0.85)

# Note: one outlier (H=1.27, n=21 transitions) excluded from clean dataset
# Shown as text note in the panel title area
ax2.text(0.5, 1.08, 'n = 15  (1 excluded)', fontsize=5.5, color='#888888',
         ha='center', va='bottom', transform=ax2.transAxes, style='italic')

# Mean and median annotations (right side)
h_vals = [r['h'] for r in clean]
ax2.annotate(f'$\\bar{{H}} = {np.mean(h_vals):.3f}$',
             xy=(0.3, np.mean(h_vals)), fontsize=6, color='#444444',
             va='center', ha='left',
             bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                       edgecolor='#CCCCCC', alpha=0.8))

ax2.set_xlim(-0.5, 1.4)
ax2.set_ylim(0.55, 1.10)
ax2.set_xticks([])
ax2.set_ylabel('$H_{\\mathrm{DFA}}$')
ax2.set_xlabel('$C_4$ recordings')

# Legend for cell types
handles_b = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_HEK,
           markeredgecolor='black', markersize=5, label='HEK293'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_MOUSE,
           markeredgecolor='black', markersize=5, label='Mouse CM'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_RABBIT,
           markeredgecolor='black', markersize=5, label='Rabbit CM'),
]
ax2.legend(handles=handles_b, loc='upper left', framealpha=0.9, fontsize=6,
           edgecolor='#CCCCCC')

ax2.text(-0.08, 1.14, 'b', transform=ax2.transAxes, fontsize=12,
         fontweight='bold', va='top')

# ═══════════════════════════════════════════════════════════════════════════
# Panel (c): H vs recording duration
# ═══════════════════════════════════════════════════════════════════════════

# Reference lines
ax3.axhline(0.833, color=COL_PRED, linestyle='--', linewidth=1.0, zorder=1)
ax3.axhline(0.69, color=COL_CEIL, linestyle=':', linewidth=0.8, zorder=1)

# Aging crossover time (from paper: t_cross ~ 47 s for C4)
ax3.axvline(47, color='#BBBBBB', linestyle='--', linewidth=0.6, zorder=0)
ax3.text(55, 0.78, '$t_{\\mathrm{cross}}$\n$\\approx 47$ s', fontsize=5.5,
         color='#999999', va='top', ha='left')

for r in clean:
    if r['cell_type'] == 'HEK293':
        color = COL_HEK
    elif r['cell_type'] == 'Mouse CM':
        color = COL_MOUSE
    else:
        color = COL_RABBIT

    msize = 5
    edgewidth = 0.5
    if r['dev'] < 2.0:
        edgewidth = 1.5
        msize = 7

    ax3.plot(r['duration'], r['h'], 'o', color=color, markersize=msize,
             markeredgewidth=edgewidth, markeredgecolor='black', zorder=10, alpha=0.85)

# Outlier
for r in outlier:
    ax3.plot(r['duration'], r['h'], 'x', color='#AAAAAA', markersize=6,
             markeredgewidth=1.0, zorder=8)

# Labels for best matches
best1 = next(r for r in clean if 'Cell01_60mV' in r['name'])
best2 = next(r for r in clean if 'Cell08_60mV' in r['name'])
ax3.annotate(f'0.4%', xy=(best1['duration'], best1['h']),
             xytext=(best1['duration'] + 40, best1['h'] + 0.04),
             fontsize=6, color='#444444',
             arrowprops=dict(arrowstyle='->', color='#888888', lw=0.6))
ax3.annotate(f'1.6%', xy=(best2['duration'], best2['h']),
             xytext=(best2['duration'] - 80, best2['h'] + 0.06),
             fontsize=6, color='#444444',
             arrowprops=dict(arrowstyle='->', color='#888888', lw=0.6))

ax3.set_xlabel('Recording duration (s)')
ax3.set_ylabel('$H_{\\mathrm{DFA}}$')
ax3.set_ylim(0.55, 1.10)
ax3.set_xlim(0, 1100)

# Text for prediction
ax3.text(600, 0.84, '$H_{\\mathrm{pred}} = 0.833$', fontsize=6,
         color=COL_PRED, ha='center', va='bottom')

ax3.text(-0.15, 1.14, 'c', transform=ax3.transAxes, fontsize=12,
         fontweight='bold', va='top')

# ── Save ───────────────────────────────────────────────────────────────────
outpath = os.path.join(os.path.dirname(__file__), 'figure4_own_dfa_results.png')
fig.savefig(outpath, dpi=300)
fig.savefig(outpath.replace('.png', '.pdf'))
print(f"Saved: {outpath}")
print(f"Saved: {outpath.replace('.png', '.pdf')}")
plt.close()
