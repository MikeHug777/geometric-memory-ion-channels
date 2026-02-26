#!/usr/bin/env python3
"""
Figure 1: The Burnside derivation in four steps
Paper: "Pore Symmetry Determines Fractal Memory in Ion Channels"

Horizontal flow diagram: (a) Binary H-bond ring → (b) Equipartition → (c) Bouchaud → (d) Lowen-Teich → H
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, Circle, FancyBboxPatch
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

COL_C4 = '#D73027'

# ── FIGURE ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 4, figsize=(7.2, 3.0))
for ax in axes:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')

# ═══════════════════════════════════════════════════════════════════════════
# Panel (a): Binary hydrogen-bond ring under C4 rotation
# ═══════════════════════════════════════════════════════════════════════════
ax = axes[0]
ax.set_title('$C_4$ H-bond ring', fontsize=8, fontweight='bold', pad=8)

# Draw ring with 4 sites
center = (0.5, 0.55)
radius = 0.25
angles = [90, 0, 270, 180]  # Top, Right, Bottom, Left
site_colors = ['black', 'white', 'black', 'white']  # Example: 1,0,1,0
labels = ['1', '0', '1', '0']

# Draw circle (ring backbone)
circle = plt.Circle(center, radius, fill=False, edgecolor='#BBBBBB',
                     linewidth=1.5, linestyle='--')
ax.add_patch(circle)

# Draw sites
for i, (angle, color, label) in enumerate(zip(angles, site_colors, labels)):
    x = center[0] + radius * np.cos(np.radians(angle))
    y = center[1] + radius * np.sin(np.radians(angle))
    site = plt.Circle((x, y), 0.06, facecolor=color, edgecolor='black',
                       linewidth=1.0, zorder=10)
    ax.add_patch(site)
    txt_color = 'white' if color == 'black' else 'black'
    ax.text(x, y, label, fontsize=8, ha='center', va='center',
            color=txt_color, fontweight='bold', zorder=11)

# Rotation arrow
from matplotlib.patches import Arc
arc = Arc(center, 0.18, 0.18, angle=0, theta1=30, theta2=330,
          color='#888888', linewidth=1.0, linestyle='-')
ax.add_patch(arc)
ax.annotate('', xy=(center[0] + 0.07, center[1] + 0.06),
            xytext=(center[0] + 0.09, center[1] + 0.03),
            arrowprops=dict(arrowstyle='->', color='#888888', lw=1.0))

# Formula
ax.text(0.5, 0.05, '$2^4 = 16$ states\n$\\downarrow$ $C_4$ rotation\n$B(4) = 6$ orbits',
        fontsize=7, ha='center', va='center', color='#444444',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#F5F5F5',
                  edgecolor='#CCCCCC'))

ax.text(0.02, 0.97, 'a', fontsize=11, fontweight='bold', va='top',
        transform=ax.transAxes)

# ═══════════════════════════════════════════════════════════════════════════
# Panel (b): Equipartition → trap depth
# ═══════════════════════════════════════════════════════════════════════════
ax = axes[1]
ax.set_title('Equipartition', fontsize=8, fontweight='bold', pad=8)

# Draw energy levels (orbits as degrees of freedom)
n_orbits = 6
for i in range(n_orbits):
    y = 0.3 + i * 0.1
    ax.plot([0.25, 0.75], [y, y], '-', color=COL_C4, linewidth=1.5, alpha=0.7)
    ax.text(0.78, y, f'orbit {i+1}', fontsize=5, va='center', color='#888888')

# Bracket
ax.annotate('', xy=(0.18, 0.3), xytext=(0.18, 0.8),
            arrowprops=dict(arrowstyle='|-|', color='black', lw=1.0))
ax.text(0.10, 0.55, '$B(n)$\ndof', fontsize=6.5, ha='center', va='center',
        rotation=90, color='#444444')

# Formula
ax.text(0.5, 0.05, '$E_0 = B(n) \\cdot kT/2$\n$kT$ cancels later',
        fontsize=7, ha='center', va='center', color='#444444',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#F5F5F5',
                  edgecolor='#CCCCCC'))

ax.text(0.02, 0.97, 'b', fontsize=11, fontweight='bold', va='top',
        transform=ax.transAxes)

# ═══════════════════════════════════════════════════════════════════════════
# Panel (c): Bouchaud trap model → power-law dwell times
# ═══════════════════════════════════════════════════════════════════════════
ax = axes[2]
ax.set_title('Bouchaud traps', fontsize=8, fontweight='bold', pad=8)

# Draw power-law distribution
tau = np.linspace(0.05, 0.95, 200)
psi = 0.25 * tau**(-1.33)  # ~tau^{-(1+mu)} with mu=1/3
psi = np.clip(psi, 0, 0.85)
ax.fill_between(tau, 0.25, psi, alpha=0.15, color=COL_C4, zorder=2)
ax.plot(tau, psi, '-', color=COL_C4, linewidth=1.5, zorder=3)

# Axes
ax.plot([0.05, 0.95], [0.25, 0.25], '-', color='black', linewidth=0.8)
ax.plot([0.05, 0.05], [0.25, 0.90], '-', color='black', linewidth=0.8)
ax.text(0.5, 0.18, '$\\tau$', fontsize=8, ha='center')
ax.text(0.00, 0.60, '$\\psi(\\tau)$', fontsize=8, ha='right', rotation=90)

# Formula
ax.text(0.5, 0.05, '$\\psi(\\tau) \\sim \\tau^{-(1+\\mu)}$\n$\\mu = 2/B(n)$',
        fontsize=7, ha='center', va='center', color='#444444',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#F5F5F5',
                  edgecolor='#CCCCCC'))

ax.text(0.08, 0.98, 'c', fontsize=11, fontweight='bold', va='top',
        transform=ax.transAxes)

# ═══════════════════════════════════════════════════════════════════════════
# Panel (d): Lowen-Teich → H
# ═══════════════════════════════════════════════════════════════════════════
ax = axes[3]
ax.set_title('Lowen-Teich', fontsize=8, fontweight='bold', pad=8)

# Fractal point process visualization
np.random.seed(42)
# Generate clustered points (fractal-like)
n_events = 40
times = np.sort(np.random.beta(0.5, 0.5, n_events))  # Clustered
for i, t in enumerate(times):
    ax.plot([t * 0.85 + 0.05, t * 0.85 + 0.05], [0.55, 0.65], '-',
            color=COL_C4, linewidth=0.8, alpha=0.7)
ax.plot([0.05, 0.90], [0.55, 0.55], '-', color='black', linewidth=0.8)
ax.text(0.50, 0.47, 'time $\\rightarrow$', fontsize=6, ha='center',
        color='#888888')
ax.text(0.50, 0.72, 'fractal gating events', fontsize=5.5, ha='center',
        color='#888888', style='italic')

# Result box with the key formula
result_box = FancyBboxPatch((0.08, 0.25), 0.84, 0.20,
                             boxstyle='round,pad=0.05',
                             facecolor='#FFF3E0', edgecolor=COL_C4,
                             linewidth=1.5, zorder=5)
ax.add_patch(result_box)
ax.text(0.50, 0.35, '$H = 1 - 1/B(n)$', fontsize=10, ha='center',
        va='center', fontweight='bold', color=COL_C4, zorder=6)

# "Parameter-free" note
ax.text(0.50, 0.12, 'Parameter-free\n(no fit parameters)', fontsize=6,
        ha='center', va='center', color='#666666', style='italic')

ax.text(0.02, 0.97, 'd', fontsize=11, fontweight='bold', va='top',
        transform=ax.transAxes)

# ── Add connecting arrows between panels ──────────────────────────────────
fig.canvas.draw()  # Force layout calculation
for i in range(3):
    pos_r = axes[i].get_position()
    pos_l = axes[i + 1].get_position()
    x_mid = (pos_r.x1 + pos_l.x0) / 2
    y_mid = (pos_r.y0 + pos_r.y1) / 2
    fig.patches.append(FancyArrowPatch(
        posA=(x_mid - 0.012, y_mid),
        posB=(x_mid + 0.012, y_mid),
        transform=fig.transFigure,
        arrowstyle='->', mutation_scale=15,
        color='#888888', linewidth=1.5, zorder=100
    ))

# Remove spines
for ax in axes:
    for spine in ax.spines.values():
        spine.set_visible(False)

plt.subplots_adjust(wspace=0.08)

# ── Save ───────────────────────────────────────────────────────────────────
outpath = os.path.join(os.path.dirname(__file__), 'figure1_burnside_derivation.png')
fig.savefig(outpath, dpi=300)
fig.savefig(outpath.replace('.png', '.pdf'))
print(f"Saved: {outpath}")
print(f"Saved: {outpath.replace('.png', '.pdf')}")
plt.close()
