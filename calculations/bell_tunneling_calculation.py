"""
Calculation SI-30a: WKB Tunneling Through H-Bond Barriers in Selectivity Filters
=================================================================================
Author:  Michael Hug
Date:    2026-02-22
Project: Interface First — Geometric Memory
Paper:   Paper 2, SI Section 30

Original contribution:
  Applies the Bell (1980) WKB tunneling correction specifically to
  the H-bond ring geometry of ion channel selectivity filters, using
  Pomes & Roux (2002) barrier parameters. Predicts that proton tunneling
  corrections Q_t are significant (>1.5) at physiological temperature for
  the 3 A pore geometry, directly linking pore structure to quantum transport.

Dependencies: numpy
"""

import numpy as np

# ============================================================
# NATURKONSTANTEN (Lehrbuch)
# ============================================================
h = 6.626e-34      # Planck-Konstante [J·s]
k_B = 1.381e-23    # Boltzmann-Konstante [J/K]
eV_to_J = 1.602e-19  # 1 eV in Joule
amu_to_kg = 1.661e-27  # 1 atomare Masseneinheit in kg

# ============================================================
# VON POMÈS & ROUX (2002) — die einzigen externen Daten
# ============================================================
# Energiebarriere zwischen zwei Wassermolekülen in Gramicidin
# Sie fanden 1-5 kcal/mol = 43-217 meV
# Wir nehmen drei Szenarien: konservativ, zentral, obere Grenze
V0_scenarios = {
    'Konservativ': 150e-3 * eV_to_J,  # 150 meV in Joule
    'Zentral':     200e-3 * eV_to_J,  # 200 meV in Joule
    'Obere Grenze':250e-3 * eV_to_J,  # 250 meV in Joule
}

# Halbe Breite der Barriere
a = 0.40e-10  # 0.40 Ångström in Meter

# ============================================================
# IONENMASSEN (Lehrbuch)
# ============================================================
ions = {
    'H⁺':  1.008 * amu_to_kg,
    'D⁺':  2.014 * amu_to_kg,
    'Li⁺': 6.941 * amu_to_kg,
    'Na⁺': 22.99 * amu_to_kg,
    'K⁺':  39.10 * amu_to_kg,
    'Cs⁺': 132.9 * amu_to_kg,
}

# ============================================================
# TEMPERATUREN
# ============================================================
temperatures = {
    '15°C': 288,   # Kelvin
    '25°C': 298,
    '37°C': 310,
    '45°C': 318,
}

print("=" * 70)
print("BELL TUNNELING MODEL — Schritt-für-Schritt-Berechnung")
print("=" * 70)

# ============================================================
# SCHRITT 0: De Broglie Wellenlängen (zum Verständnis)
# ============================================================
print("\n--- SCHRITT 0: De Broglie Wellenlängen bei 37°C (310 K) ---")
print("Formel: λ = h / √(2 m kT)  [de Broglie bei thermischer Energie E = kT]")
print()
T = 310
for name, mass in ions.items():
    if name == 'D⁺':
        continue
    lam = h / np.sqrt(2 * mass * k_B * T)
    pore_d = 4.0e-10  # 4.0 Å Gramicidin-Pore
    ratio = lam / pore_d * 100
    print(f"  {name:4s}: λ = {lam*1e10:.2f} Å  →  {ratio:.1f}% der Pore (4.0 Å)")

print()
print("  → H⁺ ist zu 44% 'so gross wie die Pore' → stark quantenmechanisch")
print("  → Cs⁺ ist zu 4% 'so gross wie die Pore' → klassisch")

# ============================================================
# SCHRITT 1-4: Bell-Tunneling für H⁺/D⁺ Isotopen-Effekt
# ============================================================
print("\n" + "=" * 70)
print("ISOTOPEN-EFFEKT: H⁺ vs. D⁺")
print("=" * 70)

for scenario_name, V0 in V0_scenarios.items():
    V0_meV = V0 / eV_to_J * 1000
    print(f"\n--- Szenario: {scenario_name} (V₀ = {V0_meV:.0f} meV, a = 0.40 Å) ---")

    for temp_name, T in temperatures.items():
        # SCHRITT 1: Imaginäre Frequenz am Barrierentop
        # ν‡ = (1/2π) × √(2V₀ / (m × a²))
        # Das ist die "Schwingungsfrequenz" die der Partikel hätte,
        # wenn er auf dem Gipfel der Barriere sitzen würde.
        # Leichtere Teilchen schwingen schneller → tunneln mehr.

        nu_H = (1/(2*np.pi)) * np.sqrt(2*V0 / (ions['H⁺'] * a**2))
        nu_D = (1/(2*np.pi)) * np.sqrt(2*V0 / (ions['D⁺'] * a**2))

        # SCHRITT 2: Tunneling-Parameter u = hν‡ / (2kT)
        # Das ist das Verhältnis von Quanten-Energie zu thermischer Energie.
        # u > 1: Tunneling wird wichtig
        # u < 1: klassisches Verhalten dominiert
        # u > π: Bell-Formel bricht zusammen (zu viel Tunneling)

        u_H = h * nu_H / (2 * k_B * T)
        u_D = h * nu_D / (2 * k_B * T)

        # SCHRITT 3: Tunneling-Korrektur Q_t = u / sin(u)
        # Q_t = 1.0 bedeutet: kein Tunneling (klassisch)
        # Q_t = 2.0 bedeutet: Teilchen kommt doppelt so schnell durch wie klassisch

        if u_H < np.pi and u_D < np.pi:
            Qt_H = u_H / np.sin(u_H)
            Qt_D = u_D / np.sin(u_D)

            # SCHRITT 4: Verhältnis
            # Die Tunneling-Korrektur ist für H⁺ grösser als für D⁺,
            # weil H⁺ leichter ist → höheres u → mehr Tunneling
            tunneling_ratio = Qt_H / Qt_D

            # Klassischer Faktor (~1.4) kommt von:
            # - D₂O ist 23% viskoser als H₂O
            # - D⁺ diffundiert langsamer (Masse-Effekt)
            classical_ratio = 1.40

            total_ratio = classical_ratio * tunneling_ratio

            if temp_name == '37°C':
                print(f"  Schritt 1: ν‡(H⁺) = {nu_H:.3e} Hz,  ν‡(D⁺) = {nu_D:.3e} Hz")
                print(f"  Schritt 2: u(H⁺)  = {u_H:.3f},          u(D⁺)  = {u_D:.3f}")
                print(f"  Schritt 3: Q_t(H⁺) = {Qt_H:.3f},        Q_t(D⁺) = {Qt_D:.3f}")
                print(f"  Schritt 4: Tunneling-Ratio = {tunneling_ratio:.3f}")
                print(f"  Schritt 5: × klassisch (1.40) = {total_ratio:.2f}")
            else:
                print(f"  {temp_name}: u(H⁺)={u_H:.3f}, Q_t(H⁺)={Qt_H:.3f}, Q_t(D⁺)={Qt_D:.3f}, "
                      f"Tunneling-Ratio={tunneling_ratio:.3f}, Total={total_ratio:.2f}")
        else:
            print(f"  {temp_name}: u(H⁺)={u_H:.3f} — WARNUNG: u > π, Bell-Formel ungültig!")

# ============================================================
# SCHRITT 5: Massenserie bei 37°C
# ============================================================
print("\n" + "=" * 70)
print("MASSENSERIE: Alle Ionen bei 37°C (V₀ = 200 meV)")
print("=" * 70)

V0 = V0_scenarios['Zentral']
T = 310

print(f"\n{'Ion':5s} {'Masse':>8s} {'λ_dB':>8s} {'u':>8s} {'Q_t':>8s} {'Excess':>10s}")
print("-" * 50)
for name, mass in ions.items():
    if name == 'D⁺':
        continue

    lam = h / np.sqrt(2 * np.pi * mass * k_B * T)
    nu = (1/(2*np.pi)) * np.sqrt(2*V0 / (mass * a**2))
    u = h * nu / (2 * k_B * T)

    if u < np.pi:
        Qt = u / np.sin(u)
        excess = f"{(Qt-1)*100:.1f}%"
    else:
        Qt = float('inf')
        excess = ">>100%"

    print(f"{name:5s} {mass/amu_to_kg:7.1f}u {lam*1e10:7.2f}Å {u:7.3f} {Qt:7.3f}  {excess:>8s}")

print()
print("Lesebeispiel:")
print("  H⁺:  Q_t = 2.02 → Proton kommt 102% schneller durch als klassisch erwartet")
print("  Li⁺: Q_t = 1.10 → Lithium kommt 10% schneller durch")
print("  Cs⁺: Q_t = 1.00 → Cäsium verhält sich rein klassisch (kein Tunneling)")

# ============================================================
# MASSENSERIE bei allen Temperaturen
# ============================================================
print("\n" + "=" * 70)
print("MASSENSERIE: Alle Ionen bei allen Temperaturen (V₀ = 200 meV)")
print("=" * 70)

for temp_name, T in temperatures.items():
    print(f"\n  {temp_name} ({T} K):")
    print(f"  {'Ion':5s} {'u':>7s} {'Q_t':>7s}")
    print(f"  {'-'*20}")
    for name, mass in ions.items():
        if name == 'D⁺':
            continue
        nu = (1/(2*np.pi)) * np.sqrt(2*V0 / (mass * a**2))
        u = h * nu / (2 * k_B * T)
        if u < np.pi:
            Qt = u / np.sin(u)
        else:
            Qt = float('inf')
        print(f"  {name:5s} {u:7.3f} {Qt:7.3f}")

print("\n" + "=" * 70)
print("ZUSAMMENFASSUNG")
print("=" * 70)
print()
print("1. Von Pomès & Roux brauchen wir: V₀ (Barrierenhöhe) und a (Barrierenbreite)")
print("2. Die Bell-Formel Q_t = u/sin(u) gibt den Tunneling-Korrekturfaktor")
print("3. Für H⁺ bei 37°C: Q_t ≈ 2.0 → starkes Tunneling")
print("4. Für Cs⁺ bei 37°C: Q_t ≈ 1.0 → kein Tunneling")
print("5. Vorhergesagtes g(H⁺)/g(D⁺) = 1.4 × (Q_t(H⁺)/Q_t(D⁺)) ≈ 1.8-2.3")
print("6. DeCoursey hat 2.0 gemessen → passt genau!")
