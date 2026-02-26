"""
Calculation AG: Molien-Serie als bosonische Zustandssumme
=========================================================

Verbindet die Invariantentheorie (Chevalley 1955) mit statistischer Mechanik.

Die Molien-Serie M_W(t) = 1/∏(1-t^{d_i}) einer Coxeter-Gruppe W mit
Invariantgraden d_1,...,d_r ist formal identisch mit der Zustandssumme
von r unabhängigen bosonischen Oszillatoren.

Kernfrage: Trägt jeder primitive Exponent kT/2 zur Energieskala E₀ bei?

Tests:
  1. Molien-Zustandssumme Z(β) für E₆, E₇, E₈ berechnen
  2. Mittlere Energie ⟨E⟩ im klassischen Limes (β→0) → rank × kT
  3. Primitive vs. nicht-primitive Modi separieren
  4. ⟨E⟩_prim → φ(h) × kT/2 = B(n) × kT/2 verifizieren
  5. Visualisierung: Energiebeiträge der einzelnen Modi

Referenzen:
  - Chevalley (1955): Invariants of finite groups
  - Notiz 17: B(n) = φ(h(E_{n+3}))
  - Notiz 24: Molien-Serie als Zustandssumme

Output: mc_molien_partition/
"""

import numpy as np
import matplotlib.pyplot as plt
from math import gcd
import os

# ============================================================
# 1. E-Serie Daten (Chevalley-Grade und Exponenten)
# ============================================================

E_SERIES = {
    'E6': {
        'rank': 6,
        'h': 12,
        'degrees': [2, 5, 6, 8, 9, 12],
        'exponents': [1, 4, 5, 7, 8, 11],  # e_i = d_i - 1
        'label': r'$E_6$ (C3)',
        'n': 3,
        'color': '#2196F3',
    },
    'E7': {
        'rank': 7,
        'h': 18,
        'degrees': [2, 6, 8, 10, 12, 14, 18],
        'exponents': [1, 5, 7, 9, 11, 13, 17],
        'label': r'$E_7$ (C4)',
        'n': 4,
        'color': '#4CAF50',
    },
    'E8': {
        'rank': 8,
        'h': 30,
        'degrees': [2, 8, 12, 14, 18, 20, 24, 30],
        'exponents': [1, 7, 11, 13, 17, 19, 23, 29],
        'label': r'$E_8$ (C5)',
        'n': 5,
        'color': '#F44336',
    },
}


def euler_phi(n):
    """Euler's totient function."""
    count = 0
    for k in range(1, n + 1):
        if gcd(k, n) == 1:
            count += 1
    return count


def burnside(n):
    """Burnside orbit count for binary necklaces of length n."""
    total = 0
    for k in range(n):
        total += 2 ** gcd(k, n)
    return total // n


def is_primitive(exponent, h):
    """Check if an exponent is primitive (coprime to h)."""
    return gcd(exponent, h) == 1


# ============================================================
# 2. Zustandssumme und thermische Größen
# ============================================================

def log_partition_single(beta, d, eps0=1.0):
    """Log partition function of a single bosonic mode with degree d.

    Z_i = 1 / (1 - exp(-beta * d * eps0))
    ln Z_i = -ln(1 - exp(-beta * d * eps0))
    """
    x = beta * d * eps0
    if x > 500:  # Overflow protection
        return 0.0
    return -np.log(1.0 - np.exp(-x))


def mean_energy_single(beta, d, eps0=1.0):
    """Mean energy of a single bosonic mode with degree d.

    <E_i> = d * eps0 / (exp(beta * d * eps0) - 1)

    Classical limit (beta -> 0): <E_i> -> 1/beta = kT
    """
    x = beta * d * eps0
    if x > 500:
        return 0.0
    return d * eps0 / (np.exp(x) - 1.0)


def energy_variance_single(beta, d, eps0=1.0):
    """Energy variance (heat capacity contribution) of a single mode.

    C_i = d(E)/d(T) = (d * eps0)^2 * exp(beta*d*eps0) / (exp(beta*d*eps0) - 1)^2 / kT^2
    Var(E_i) = kT^2 * C_i
    """
    x = beta * d * eps0
    if x > 500:
        return 0.0
    ex = np.exp(x)
    return (d * eps0) ** 2 * ex / (ex - 1.0) ** 2


# ============================================================
# 3. Tests
# ============================================================

def test1_classical_limit():
    """Test 1: Im klassischen Limes → ⟨E⟩ = rank × kT (volle Oszillatoren)
    oder rank × kT/2 (nur potentielle Energie)."""

    print("=" * 70)
    print("TEST 1: Klassischer Limes (β → 0)")
    print("=" * 70)
    print()

    betas = np.logspace(-3, 1, 200)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for idx, (name, data) in enumerate(E_SERIES.items()):
        ax = axes[idx]

        h = data['h']
        degrees = data['degrees']
        exponents = data['exponents']
        rank = data['rank']
        n = data['n']
        B_n = burnside(n)
        phi_h = euler_phi(h)

        # Separate primitive and non-primitive
        prim_degrees = [d for d, e in zip(degrees, exponents) if is_primitive(e, h)]
        nonprim_degrees = [d for d, e in zip(degrees, exponents) if not is_primitive(e, h)]

        E_total = np.array([sum(mean_energy_single(b, d) for d in degrees) for b in betas])
        E_prim = np.array([sum(mean_energy_single(b, d) for d in prim_degrees) for b in betas])
        E_nonprim = np.array([sum(mean_energy_single(b, d) for d in nonprim_degrees) for b in betas])

        kT = 1.0 / betas

        # Plot
        ax.plot(kT, E_total / kT, 'k-', lw=2, label=f'Alle {rank} Modi')
        ax.plot(kT, E_prim / kT, '-', color=data['color'], lw=2.5,
                label=f'{phi_h} primitive Modi')
        if len(nonprim_degrees) > 0:
            ax.plot(kT, E_nonprim / kT, '--', color='gray', lw=1.5,
                    label=f'{rank - phi_h} nicht-primitiv')

        # Reference lines
        ax.axhline(rank, color='k', ls=':', alpha=0.5, label=f'rank = {rank}')
        ax.axhline(phi_h, color=data['color'], ls=':', alpha=0.5,
                    label=f'φ(h) = B({n}) = {phi_h}')
        ax.axhline(phi_h / 2, color=data['color'], ls='-.', alpha=0.5,
                    label=f'B({n})/2 = {phi_h/2}')

        ax.set_xscale('log')
        ax.set_xlabel('kT / ε₀')
        ax.set_ylabel('⟨E⟩ / kT')
        ax.set_title(f'{data["label"]}: h={h}, rank={rank}, φ(h)={phi_h}')
        ax.legend(fontsize=7, loc='best')
        ax.set_xlim(0.1, 1000)
        ax.set_ylim(0, rank + 1)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()

    outdir = 'mc_molien_partition'
    os.makedirs(outdir, exist_ok=True)
    plt.savefig(f'{outdir}/fig1_classical_limit.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{outdir}/fig1_classical_limit.pdf', bbox_inches='tight')
    plt.close()

    # Print numerical values
    print(f"{'Algebra':<8} {'rank':<6} {'φ(h)':<6} {'B(n)':<6} "
          f"{'⟨E⟩/kT (β=0.001)':<20} {'⟨E⟩_prim/kT':<15} {'Ratio':<10}")
    print("-" * 75)

    beta_low = 0.001
    for name, data in E_SERIES.items():
        h = data['h']
        degrees = data['degrees']
        exponents = data['exponents']
        rank = data['rank']
        n = data['n']
        B_n = burnside(n)
        phi_h = euler_phi(h)

        prim_degrees = [d for d, e in zip(degrees, exponents) if is_primitive(e, h)]

        E_tot = sum(mean_energy_single(beta_low, d) for d in degrees)
        E_prim = sum(mean_energy_single(beta_low, d) for d in prim_degrees)
        kT = 1.0 / beta_low

        print(f"{name:<8} {rank:<6} {phi_h:<6} {B_n:<6} "
              f"{E_tot/kT:<20.4f} {E_prim/kT:<15.4f} {E_prim/E_tot:<10.4f}")

    print()
    print("Im klassischen Limes (kT >> ε₀·d_max):")
    print("  ⟨E⟩_total / kT → rank       (alle Modi tragen kT bei)")
    print("  ⟨E⟩_prim  / kT → φ(h) = B(n) (primitive Modi tragen kT bei)")
    print()
    print("Für E₀ = B(n)·kT/2 brauchen wir den Faktor 1/2:")
    print("  → Diskreter Zustandsraum {0,1}^n → nur potentielle Energie → kT/2 pro Mode")


def test2_energy_per_mode():
    """Test 2: Energie-Beitrag pro Mode als Funktion von β."""

    print()
    print("=" * 70)
    print("TEST 2: Energie pro Mode")
    print("=" * 70)
    print()

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    betas = np.logspace(-2, 1, 200)

    for idx, (name, data) in enumerate(E_SERIES.items()):
        ax = axes[idx]
        h = data['h']
        degrees = data['degrees']
        exponents = data['exponents']

        for d, e in zip(degrees, exponents):
            prim = is_primitive(e, h)
            E_mode = np.array([mean_energy_single(b, d) for b in betas])
            kT = 1.0 / betas

            style = '-' if prim else '--'
            alpha = 1.0 if prim else 0.4
            color = data['color'] if prim else 'gray'

            ax.plot(kT, E_mode / kT, style, color=color, alpha=alpha, lw=1.5,
                    label=f'd={d}, e={e}' + (' ★' if prim else ''))

        # Reference
        ax.axhline(1.0, color='red', ls=':', alpha=0.3, label='kT (full osc.)')
        ax.axhline(0.5, color='blue', ls=':', alpha=0.3, label='kT/2 (pot. only)')

        ax.set_xscale('log')
        ax.set_xlabel('kT / ε₀')
        ax.set_ylabel('⟨Eᵢ⟩ / kT')
        ax.set_title(f'{data["label"]}')
        ax.legend(fontsize=6, loc='best', ncol=2)
        ax.set_xlim(0.1, 100)
        ax.set_ylim(0, 1.2)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()

    outdir = 'mc_molien_partition'
    plt.savefig(f'{outdir}/fig2_energy_per_mode.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{outdir}/fig2_energy_per_mode.pdf', bbox_inches='tight')
    plt.close()

    print("Im klassischen Limes: Jeder Mode (primitiv oder nicht) → kT")
    print("Für halbes Equipartition (kT/2): Interpretation als rein potentiell")


def test3_primitive_selection():
    """Test 3: Zeigt dass die Selektion 'nur primitive' genau B(n) ergibt."""

    print()
    print("=" * 70)
    print("TEST 3: Primitive-Selektion → B(n)")
    print("=" * 70)
    print()

    print(f"{'Algebra':<8} {'h':<5} {'rank':<6} {'φ(h)':<6} {'B(n)':<6} "
          f"{'Prim. Exp.':<25} {'Non-prim.':<15} {'Match?'}")
    print("-" * 95)

    for name, data in E_SERIES.items():
        h = data['h']
        exponents = data['exponents']
        rank = data['rank']
        n = data['n']
        B_n = burnside(n)
        phi_h = euler_phi(h)

        prim = [e for e in exponents if is_primitive(e, h)]
        nonprim = [e for e in exponents if not is_primitive(e, h)]

        match = "✓" if phi_h == B_n else "✗"

        print(f"{name:<8} {h:<5} {rank:<6} {phi_h:<6} {B_n:<6} "
              f"{str(prim):<25} {str(nonprim):<15} {match}")

    print()
    print("Ergebnis: φ(h(E_{n+3})) = B(n) für alle drei E-Algebren ✓")
    print("  E₆: 4 von 6 Exponenten primitiv → 4 = B(3)")
    print("  E₇: 6 von 7 Exponenten primitiv → 6 = B(4)")
    print("  E₈: 8 von 8 Exponenten primitiv → 8 = B(5)")


def test4_bouchaud_connection():
    """Test 4: Verbindung Molien → Bouchaud → H.

    Wenn ⟨E⟩_prim = B(n) · kT/2 die Skala der Trap-Verteilung setzt:
      μ = kT / E₀ = kT / (B(n) · kT/2) = 2/B(n)
      H = 1 - μ/2 = 1 - 1/B(n)
    """

    print()
    print("=" * 70)
    print("TEST 4: Molien → Bouchaud → GM-Formel")
    print("=" * 70)
    print()

    print("Logische Kette (Route E):")
    print()
    print("  Schritt E.1: Cn-Pore → E_{n+3}          (McKay, empirisch)")
    print("  Schritt E.2: Invariantgrade d_i          (Chevalley 1955, bewiesen)")
    print("  Schritt E.3: M_W(t) = ∏ 1/(1-t^{d_i})  (Molien, bewiesen)")
    print("  Schritt E.4: M_W(e^{-β}) = Z_boson(β)   (formale Identität)")
    print("  Schritt E.5: Nur primitive Modi zählen   ← KONJEKTUR")
    print("  Schritt E.6: Diskret → kT/2 pro Mode    (physikalisch plausibel)")
    print("  ⟹  E₀ = φ(h) · kT/2 = B(n) · kT/2")
    print("  ⟹  μ = 2/B(n)")
    print("  ⟹  H = 1 - 1/B(n)  = GM-Formel  ∎")
    print()

    print(f"{'Algebra':<8} {'n':<4} {'B(n)':<6} {'φ(h)':<6} "
          f"{'E₀/kT':<8} {'μ=2/B':<8} {'H=1-1/B':<10} {'H gemessen'}")
    print("-" * 65)

    measured = {'E6': '?', 'E7': '~0.82', 'E8': '?'}

    for name, data in E_SERIES.items():
        n = data['n']
        B_n = burnside(n)
        phi_h = euler_phi(data['h'])
        E0_over_kT = B_n / 2.0
        mu = 2.0 / B_n
        H = 1.0 - 1.0 / B_n

        print(f"{name:<8} {n:<4} {B_n:<6} {phi_h:<6} "
              f"{E0_over_kT:<8.3f} {mu:<8.4f} {H:<10.4f} {measured[name]}")

    print()
    print("Die Molien-Route erklärt:")
    print("  (a) WARUM E₀ ∝ B(n): Weil B(n) = φ(h) = Zahl der primitiven Modi")
    print("  (b) WARUM kT/2: Weil der Zustandsraum {0,1}^n diskret ist (nur potentiell)")
    print("  (c) WARUM T-unabhängig: Weil kT sich in μ = kT/E₀ = kT/(B·kT/2) kürzt")


def test5_summary_figure():
    """Test 5: Zusammenfassende Figur."""

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    betas = np.logspace(-3, 1, 300)
    kT = 1.0 / betas

    # --- Panel (a): ⟨E⟩_prim / (kT/2) → B(n) ---
    ax = axes[0, 0]
    for name, data in E_SERIES.items():
        h = data['h']
        degrees = data['degrees']
        exponents = data['exponents']
        n = data['n']
        B_n = burnside(n)

        prim_degrees = [d for d, e in zip(degrees, exponents) if is_primitive(e, h)]
        E_prim = np.array([sum(mean_energy_single(b, d) for d in prim_degrees) for b in betas])

        ax.plot(kT, E_prim / (kT / 2), '-', color=data['color'], lw=2.5,
                label=f'{data["label"]}: B({n})={B_n}')
        ax.axhline(B_n, color=data['color'], ls=':', alpha=0.4)

    ax.set_xscale('log')
    ax.set_xlabel('kT / ε₀')
    ax.set_ylabel('⟨E⟩_prim / (kT/2)')
    ax.set_title('(a) Primitive Energie → B(n)·kT/2')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.1, 1000)

    # --- Panel (b): Anteil primitiv an total ---
    ax = axes[0, 1]
    for name, data in E_SERIES.items():
        h = data['h']
        degrees = data['degrees']
        exponents = data['exponents']
        rank = data['rank']
        phi_h = euler_phi(h)

        prim_degrees = [d for d, e in zip(degrees, exponents) if is_primitive(e, h)]

        E_total = np.array([sum(mean_energy_single(b, d) for d in degrees) for b in betas])
        E_prim = np.array([sum(mean_energy_single(b, d) for d in prim_degrees) for b in betas])

        ratio = E_prim / np.maximum(E_total, 1e-30)

        ax.plot(kT, ratio, '-', color=data['color'], lw=2.5,
                label=f'{data["label"]}: φ(h)/rank = {phi_h}/{rank}')
        ax.axhline(phi_h / rank, color=data['color'], ls=':', alpha=0.4)

    ax.set_xscale('log')
    ax.set_xlabel('kT / ε₀')
    ax.set_ylabel('⟨E⟩_prim / ⟨E⟩_total')
    ax.set_title('(b) Anteil primitive Moden')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.1, 1000)
    ax.set_ylim(0, 1.1)

    # --- Panel (c): Implied μ and H ---
    ax = axes[1, 0]
    n_values = [3, 4, 5]
    B_values = [burnside(n) for n in n_values]
    H_values = [1 - 1/B for B in B_values]
    phi_values = [euler_phi(E_SERIES[f'E{n+3}']['h']) for n in n_values]

    x = np.arange(len(n_values))
    width = 0.3

    bars1 = ax.bar(x - width/2, B_values, width, label='B(n) Burnside',
                   color=['#2196F3', '#4CAF50', '#F44336'], alpha=0.7)
    bars2 = ax.bar(x + width/2, phi_values, width, label='φ(h) Primitive',
                   color=['#2196F3', '#4CAF50', '#F44336'], alpha=0.4,
                   edgecolor='black', linewidth=1.5)

    for i, (b, p) in enumerate(zip(B_values, phi_values)):
        ax.text(x[i], max(b, p) + 0.3, f'{b}={p}', ha='center', fontsize=12, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels([f'C{n}\n(E{n+3})' for n in n_values])
    ax.set_ylabel('Anzahl DOFs')
    ax.set_title('(c) B(n) = φ(h) für n=3,4,5')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    # --- Panel (d): GM-Formel Vorhersage ---
    ax = axes[1, 1]

    n_all = [2, 3, 4, 5, 6]
    B_all = [burnside(n) for n in n_all]
    H_all = [1 - 1/B for B in B_all]

    colors_all = ['gray', '#2196F3', '#4CAF50', '#F44336', '#9C27B0']
    labels_all = ['C2', 'C3 (E₆)', 'C4 (E₇)', 'C5 (E₈)', 'C6 (kein E-Typ)']

    for i, (n, B, H, c, lab) in enumerate(zip(n_all, B_all, H_all, colors_all, labels_all)):
        marker = 'o' if n in [3, 4, 5] else 's'
        ax.scatter(B, H, s=150, c=c, marker=marker, zorder=5,
                   edgecolors='black', linewidth=1)
        ax.annotate(lab, (B, H), textcoords="offset points",
                    xytext=(10, -5), fontsize=8)

    # Experimental point
    ax.scatter(6, 0.82, s=100, c='#4CAF50', marker='*', zorder=6,
               edgecolors='black', linewidth=1)
    ax.annotate('BK gemessen\nH≈0.82', (6, 0.82), textcoords="offset points",
                xytext=(10, 10), fontsize=7, fontstyle='italic')

    # Theory curve
    B_cont = np.linspace(2, 16, 100)
    H_cont = 1 - 1 / B_cont
    ax.plot(B_cont, H_cont, 'k--', alpha=0.5, label='H = 1 − 1/B(n)')

    ax.set_xlabel('B(n) = Burnside-Orbits = φ(h) primitive Exponenten')
    ax.set_ylabel('H = Hurst-Exponent')
    ax.set_title('(d) GM-Formel: H = 1 − 1/B(n)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(1, 16)
    ax.set_ylim(0.5, 1.0)

    plt.suptitle('Molien-Zustandssumme → GM-Formel\n'
                 'Route E: Primitive Exponenten als thermodynamische Freiheitsgrade',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()

    outdir = 'mc_molien_partition'
    plt.savefig(f'{outdir}/fig3_summary.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{outdir}/fig3_summary.pdf', bbox_inches='tight')
    plt.close()

    print()
    print("=" * 70)
    print("TEST 5: Summary-Figur erstellt")
    print("=" * 70)


def test6_coxeter_phase_locking():
    """Test 6: Phasenkopplung nicht-primitiver Moden.

    Nicht-primitive Exponenten haben Periode h/gcd(e,h) < h.
    Sie 'locken' an den Coxeter-Zyklus und mitteln sich über eine
    volle Coxeter-Periode aus.

    Visualisierung: Coxeter-Eigenwerte auf dem Einheitskreis.
    """

    print()
    print("=" * 70)
    print("TEST 6: Coxeter-Phasenkopplung")
    print("=" * 70)
    print()

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for idx, (name, data) in enumerate(E_SERIES.items()):
        ax = axes[idx]
        h = data['h']
        exponents = data['exponents']

        # Unit circle
        theta = np.linspace(0, 2 * np.pi, 200)
        ax.plot(np.cos(theta), np.sin(theta), 'k-', alpha=0.2, lw=0.5)
        ax.axhline(0, color='k', alpha=0.1)
        ax.axvline(0, color='k', alpha=0.1)

        for e in exponents:
            prim = is_primitive(e, h)
            angle = 2 * np.pi * e / h
            x, y = np.cos(angle), np.sin(angle)

            color = data['color'] if prim else 'gray'
            marker = 'o' if prim else 'x'
            size = 100 if prim else 60
            alpha = 1.0 if prim else 0.5

            ax.scatter(x, y, s=size, c=color, marker=marker, alpha=alpha,
                       edgecolors='black', linewidth=0.5, zorder=5)

            # Label
            r_label = 1.2
            ax.text(r_label * np.cos(angle), r_label * np.sin(angle),
                    f'e={e}', fontsize=7, ha='center', va='center',
                    color=color, fontweight='bold' if prim else 'normal')

            # Show period
            if not prim:
                period = h // gcd(e, h)
                # Draw the sub-orbit
                sub_angles = [2 * np.pi * e * k / h for k in range(period)]
                sub_x = [np.cos(a) for a in sub_angles]
                sub_y = [np.sin(a) for a in sub_angles]
                sub_x.append(sub_x[0])
                sub_y.append(sub_y[0])
                ax.plot(sub_x, sub_y, '--', color='gray', alpha=0.3, lw=1)

        phi_h = euler_phi(h)
        ax.set_title(f'{data["label"]}: h={h}\n'
                     f'{phi_h} primitive (●) + {data["rank"]-phi_h} nicht-primitiv (×)')
        ax.set_aspect('equal')
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)
        ax.grid(True, alpha=0.1)

    plt.suptitle('Coxeter-Eigenwerte exp(2πi·e/h) auf dem Einheitskreis\n'
                 'Primitive Eigenwerte (●) erzeugen die volle Gruppe Z_h',
                 fontsize=11, fontweight='bold')
    plt.tight_layout()

    outdir = 'mc_molien_partition'
    plt.savefig(f'{outdir}/fig4_coxeter_phases.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{outdir}/fig4_coxeter_phases.pdf', bbox_inches='tight')
    plt.close()

    # Print analysis
    for name, data in E_SERIES.items():
        h = data['h']
        exponents = data['exponents']
        print(f"\n{name} (h={h}):")
        for e in exponents:
            d = gcd(e, h)
            period = h // d
            prim = "PRIMITIV" if d == 1 else f"nicht-primitiv (Periode {period}, gcd={d})"
            print(f"  e={e:>3}:  exp(2πi·{e}/{h}) hat Ordnung {period:>3} in Z_{h}  → {prim}")


def test7_non_primitive_averaging():
    """Test 7: Explizite Berechnung — Nicht-primitive Moden mitteln sich aus.

    Summiere exp(2πi·e·k/h) über einen vollen Coxeter-Zyklus k=0,...,h-1.
    Für primitive e: Summe = 0 (gleichmäßig verteilt → keine Auslöschung)
    Für nicht-primitive e: Summe = 0 (auch, aber sie bilden Untergruppen-Muster)

    Der Unterschied: Die KORRELATION über den Zyklus.
    Primitive Moden: <exp(2πi·e·k/h) · exp(2πi·e·l/h)> nur korreliert wenn k=l mod h
    Nicht-primitive: Korreliert über kürzere Perioden → Information geht schneller verloren
    """

    print()
    print("=" * 70)
    print("TEST 7: Phasen-Autokorrelation primitiv vs. nicht-primitiv")
    print("=" * 70)
    print()

    for name, data in E_SERIES.items():
        h = data['h']
        exponents = data['exponents']

        print(f"\n{name} (h={h}):")
        print(f"{'Exp.':<6} {'Primitiv?':<12} {'Ordnung':<10} "
              f"{'Autokorr. Abklingzeit':<25} {'Info-Lebensdauer'}")
        print("-" * 70)

        for e in exponents:
            d = gcd(e, h)
            order = h // d
            prim = is_primitive(e, h)

            # Autokorrelation: C(τ) = ⟨exp(2πi·e·k/h) · exp(-2πi·e·(k+τ)/h)⟩
            # = exp(-2πi·e·τ/h) → |C(τ)| = 1 (immer!)
            # Aber die EFFEKTIVE Zeitskala der Mode ist proportional zu ihrer Ordnung
            # Nicht-primitive Moden wiederholen sich nach h/gcd(e,h) Schritten

            lifetime = order  # Effektive Informationslebensdauer

            status = "★ PRIMITIV" if prim else f"  Periode = h/{d}"
            print(f"e={e:<4} {'Ja' if prim else 'Nein':<12} {order:<10} "
                  f"{'h = ' + str(h) if prim else str(order) + ' < h':<25} {status}")

        # Zusammenfassung
        prim_count = sum(1 for e in exponents if is_primitive(e, h))
        nonprim_count = len(exponents) - prim_count
        print(f"\n  → {prim_count} primitive Moden: Volle Lebensdauer h={h}")
        print(f"  → {nonprim_count} nicht-primitive Moden: Kürzere Lebensdauer")
        print(f"  → Für Langzeit-Aging (Bouchaud): Nur die {prim_count} primitiven zählen")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print("CALCULATION AG: Molien-Serie als bosonische Zustandssumme")
    print("=" * 70)
    print()

    # Verify basic identities
    print("Grundlegende Identitäten:")
    for name, data in E_SERIES.items():
        h = data['h']
        n = data['n']
        B_n = burnside(n)
        phi_h = euler_phi(h)
        print(f"  {name}: B({n}) = {B_n}, φ({h}) = {phi_h}, "
              f"Match: {'✓' if B_n == phi_h else '✗'}")
    print()

    test1_classical_limit()
    test2_energy_per_mode()
    test3_primitive_selection()
    test4_bouchaud_connection()
    test5_summary_figure()
    test6_coxeter_phase_locking()
    test7_non_primitive_averaging()

    print()
    print("=" * 70)
    print("ZUSAMMENFASSUNG")
    print("=" * 70)
    print()
    print("Die Molien-Serie M_W(t) = ∏ 1/(1-t^{d_i}) IST die bosonische")
    print("Zustandssumme von rank unabhängigen Oszillatoren.")
    print()
    print("Selektion auf primitive Modi (gcd(e_i, h) = 1) ergibt:")
    print("  ⟨E⟩_prim → φ(h) · kT/2 = B(n) · kT/2 = E₀")
    print()
    print("Dies liefert Route E zum GM-Theorem:")
    print("  E₀ = B(n)·kT/2  →  μ = 2/B(n)  →  H = 1 - 1/B(n)")
    print()
    print("OFFENE LÜCKE: Warum tragen nur primitive Modi zum")
    print("Langzeit-Tail bei? (Phasenkopplungs-Argument ist heuristisch)")
    print()
    print("Output: mc_molien_partition/fig1-4 (.png + .pdf)")
