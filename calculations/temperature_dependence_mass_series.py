"""
Calculation SI-30b: Mass-Dependent KIE Series Across Isotope Substitutions in Cn Pores
=======================================================================================
Author:  Michael Hug
Date:    2026-02-22
Project: Interface First — Geometric Memory
Paper:   Paper 2, SI Section 30

Original contribution:
  Predicts the full temperature x mass matrix of tunneling corrections
  for the biologically relevant ion series (H+, Li+, Na+, K+, Cs+),
  showing that the KIE divergence between light and heavy ions grows
  dramatically below 288 K. Provides quantitative predictions for
  cryogenic experiments that could detect tunneling even in heavier ions.

Dependencies: numpy
"""

import numpy as np

h = 6.626e-34
k_B = 1.381e-23
eV_to_J = 1.602e-19
amu_to_kg = 1.661e-27

V0 = 200e-3 * eV_to_J  # 200 meV (zentrales Szenario)
a = 0.40e-10  # 0.40 Å

ions = {
    'H⁺':  1.008 * amu_to_kg,
    'Li⁺': 6.941 * amu_to_kg,
    'Na⁺': 22.99 * amu_to_kg,
    'K⁺':  39.10 * amu_to_kg,
    'Cs⁺': 132.9 * amu_to_kg,
}

# ============================================================
# Teil 1: Biologisch relevante Temperaturen (15-45°C)
# ============================================================
print("=" * 70)
print("TEIL 1: Biologische Temperaturen — Q_t (Tunneling-Korrektur)")
print("=" * 70)
print()
print(f"{'Temp':>8s}", end="")
for name in ions:
    print(f"  {name:>8s}", end="")
print()
print("-" * 58)

for T in [288, 298, 310, 318]:
    label = f"{T-273}°C"
    print(f"{label:>8s}", end="")
    for name, mass in ions.items():
        nu = (1/(2*np.pi)) * np.sqrt(2*V0 / (mass * a**2))
        u = h * nu / (2 * k_B * T)
        if u < np.pi:
            Qt = u / np.sin(u)
        else:
            Qt = float('inf')
        print(f"  {Qt:8.3f}", end="")
    print()

print()
print("→ Bei 15°C ist H⁺ Q_t = 2.31 vs. bei 45°C Q_t = 1.94")
print("  Das ist ein MESSBARER Unterschied!")

# ============================================================
# Teil 2: Erweiterte Temperaturen — was passiert bei Kälte?
# ============================================================
print()
print("=" * 70)
print("TEIL 2: Erweiterte Temperaturen — was passiert bei Kälte?")
print("=" * 70)
print()
print("ACHTUNG: Unter 0°C friert Wasser! Gramicidin in Bilayer")
print("funktioniert nur in flüssiger Phase. Aber:")
print("- Mit Frostschutz (Glycerol) bis ca. -20°C möglich")
print("- Oder: nicht-wässrige Lösungsmittel")
print()

print(f"{'Temp':>8s}", end="")
for name in ions:
    print(f"  {name:>8s}", end="")
print()
print("-" * 58)

extended_temps = [
    (253, "-20°C"),
    (263, "-10°C"),
    (273, "0°C"),
    (283, "10°C"),
    (293, "20°C"),
    (310, "37°C"),
    (318, "45°C"),
]

for T, label in extended_temps:
    print(f"{label:>8s}", end="")
    for name, mass in ions.items():
        nu = (1/(2*np.pi)) * np.sqrt(2*V0 / (mass * a**2))
        u = h * nu / (2 * k_B * T)
        if u < np.pi:
            Qt = u / np.sin(u)
        else:
            Qt = float('inf')
        if Qt > 100:
            print(f"  {'>>100':>8s}", end="")
        else:
            print(f"  {Qt:8.3f}", end="")
    print()

# ============================================================
# Teil 3: Ab welcher Temperatur wird Li⁺ Tunneling messbar?
# ============================================================
print()
print("=" * 70)
print("TEIL 3: Ab wann wird Tunneling für jedes Ion >5% (Q_t > 1.05)?")
print("=" * 70)
print()

for name, mass in ions.items():
    nu = (1/(2*np.pi)) * np.sqrt(2*V0 / (mass * a**2))
    # Suche T wo Q_t = 1.05
    # u/sin(u) = 1.05 → u ≈ 0.555 (numerisch)
    # u = hν/(2kT) → T = hν/(2ku)
    # Für Q_t = 1.05, u ≈ 0.555
    # Für Q_t = 1.10, u ≈ 0.766

    # Numerische Suche
    found_5pct = None
    found_10pct = None
    for T_test in range(600, 10, -1):  # 600K bis 10K
        u_test = h * nu / (2 * k_B * T_test)
        if u_test < np.pi:
            Qt_test = u_test / np.sin(u_test)
            if Qt_test >= 1.05 and found_5pct is None:
                found_5pct = T_test
            if Qt_test >= 1.10 and found_10pct is None:
                found_10pct = T_test

    print(f"  {name:4s}: Q_t > 1.05 (5% Excess) ab T < {found_5pct}K ({found_5pct-273}°C)")
    if found_10pct:
        print(f"        Q_t > 1.10 (10% Excess) ab T < {found_10pct}K ({found_10pct-273}°C)")
    else:
        print(f"        Q_t > 1.10 (10% Excess): nicht erreichbar mit Bell-Modell")

# ============================================================
# Teil 4: Isotopen-Effekt g(H⁺)/g(D⁺) bei verschiedenen T
# ============================================================
print()
print("=" * 70)
print("TEIL 4: Isotopen-Effekt g(H⁺)/g(D⁺) bei verschiedenen Temperaturen")
print("=" * 70)
print()

m_H = 1.008 * amu_to_kg
m_D = 2.014 * amu_to_kg

print(f"{'Temp':>8s} {'u(H⁺)':>8s} {'u(D⁺)':>8s} {'Qt(H⁺)':>8s} {'Qt(D⁺)':>8s} {'Ratio':>8s} {'Total':>8s}")
print("-" * 60)

for T, label in extended_temps:
    nu_H = (1/(2*np.pi)) * np.sqrt(2*V0 / (m_H * a**2))
    nu_D = (1/(2*np.pi)) * np.sqrt(2*V0 / (m_D * a**2))
    u_H = h * nu_H / (2 * k_B * T)
    u_D = h * nu_D / (2 * k_B * T)

    if u_H < np.pi and u_D < np.pi:
        Qt_H = u_H / np.sin(u_H)
        Qt_D = u_D / np.sin(u_D)
        ratio = Qt_H / Qt_D
        total = 1.40 * ratio
        print(f"{label:>8s} {u_H:8.3f} {u_D:8.3f} {Qt_H:8.3f} {Qt_D:8.3f} {ratio:8.3f} {total:8.2f}")
    else:
        print(f"{label:>8s} {u_H:8.3f} — Bell-Modell bricht zusammen (u > π)!")

print()
print("=" * 70)
print("FAZIT")
print("=" * 70)
print()
print("1. JA — bei Kälte werden die Unterschiede grösser.")
print("   H⁺/D⁺ Ratio steigt von 2.00 (45°C) auf 2.21 (15°C)")
print()
print("2. Für Li⁺: 10% Excess schon bei normalen Temperaturen messbar!")
print("   Aber die Messung ist schwierig weil Li⁺ sehr schlecht")
print("   durch Gramicidin leitet (niedrige absolute Leitfähigkeit)")
print()
print("3. Für Na⁺/K⁺: Selbst bei -20°C bleibt der Effekt <5%.")
print("   → Für diese Ionen reicht der Bell-Tunneling-Effekt NICHT")
print("   → Man bräuchte Spektroskopie (Draft E) oder engere Poren")
print()
print("4. PROBLEM bei tiefen Temperaturen:")
print("   - Wasser friert ein (Bilayer bricht zusammen)")
print("   - Lipid-Bilayer werden steif unter ~5°C")
print("   - Gramicidin-Kanal-Dynamik ändert sich")
print("   → Biologisch relevanter Bereich: 15-45°C")
