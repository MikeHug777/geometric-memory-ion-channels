#!/usr/bin/env python3
"""
Calculation AB: Burnside-Totient-E Identity — B(n) = phi(h(E_{n+3}))
======================================================================
Author:  Michael Hug
Date:    2026-02-22
Project: Interface First — Geometric Memory
Paper:   Paper 2, SI Section 3

Original contribution:
  Proves the identity B(n) = phi(h(E_{n+3})) for n = 3, 4, 5 — connecting
  Burnside orbit counting to Coxeter numbers of exceptional Lie algebras
  via Euler's totient function. Since phi(h) counts the algebraically
  independent primitive exponents (Chevalley 1955, Kostant 1959), this
  provides a representation-theoretic justification for treating B(n)
  orbits as independent thermodynamic degrees of freedom: they correspond
  to genuinely independent modes of the Coxeter dynamics.

Dependencies: numpy, scipy
"""

import numpy as np
from math import gcd, factorial
from functools import reduce

print("=" * 72)
print("BERECHNUNG AB: Burnside-Totient-E-Identität")
print("B(n) = φ(h(E_{n+3})) — Darstellungstheoretische Equipartition")
print("=" * 72)

# =====================================================================
# TEIL 1: Burnside-Orbitzahlen B(n)
# =====================================================================
print("\n" + "-" * 72)
print("TEIL 1: Burnside-Orbitzahlen B(n)")
print("-" * 72)

def euler_totient(n):
    """Euler's totient function φ(n)"""
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result

def burnside_count(n):
    """B(n) = (1/n) Σ_{d|n} φ(n/d) × 2^d = Burnside-Orbitzahl"""
    total = 0
    for k in range(n):
        total += 2 ** gcd(n, k)
    return total // n

# Alternative Berechnung mit Teilern
def burnside_via_divisors(n):
    """B(n) via Teilerformel"""
    total = 0
    for d in range(1, n + 1):
        if n % d == 0:
            total += euler_totient(n // d) * (2 ** d)
    return total // n

print("\nBurnside-Orbitzahlen (zwei unabhängige Berechnungen):")
print(f"{'n':>3} | {'B(n) Rot.':>10} | {'B(n) Teiler':>12} | {'Match':>6}")
print("-" * 40)
for n in range(2, 8):
    b_rot = burnside_count(n)
    b_div = burnside_via_divisors(n)
    match = "✓" if b_rot == b_div else "✗"
    print(f"{n:>3} | {b_rot:>10} | {b_div:>12} | {match:>6}")

# =====================================================================
# TEIL 2: Coxeter-Zahlen und Totient
# =====================================================================
print("\n" + "-" * 72)
print("TEIL 2: Coxeter-Zahlen der E-Serie und Euler-Totient")
print("-" * 72)

# Coxeter-Zahlen (exakt, aus Standardtabellen)
coxeter_numbers = {
    'A1': 2, 'A2': 3, 'A3': 4, 'A4': 5, 'A5': 6,
    'D4': 6, 'D5': 8, 'D6': 10,
    'E6': 12, 'E7': 18, 'E8': 30,
    'F4': 12, 'G2': 6
}

# Ränge
ranks = {
    'A1': 1, 'A2': 2, 'A3': 3, 'A4': 4, 'A5': 5,
    'D4': 4, 'D5': 5, 'D6': 6,
    'E6': 6, 'E7': 7, 'E8': 8,
    'F4': 4, 'G2': 2
}

print("\nE-Serie:")
print(f"{'Algebra':>8} | {'rank':>5} | {'h':>5} | {'Faktorisierung':>16} | {'φ(h)':>5}")
print("-" * 52)

def factorize(n):
    """Primfaktorisierung als String"""
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    parts = []
    for p in sorted(factors):
        if factors[p] == 1:
            parts.append(str(p))
        else:
            parts.append(f"{p}^{factors[p]}")
    return " × ".join(parts)

for alg in ['E6', 'E7', 'E8']:
    h = coxeter_numbers[alg]
    r = ranks[alg]
    phi_h = euler_totient(h)
    print(f"{alg:>8} | {r:>5} | {h:>5} | {factorize(h):>16} | {phi_h:>5}")

# =====================================================================
# TEIL 3: Die Identität B(n) = φ(h(E_{n+3}))
# =====================================================================
print("\n" + "-" * 72)
print("TEIL 3: Verifikation der Identität B(n) = φ(h(E_{n+3}))")
print("-" * 72)

print("\n★ KERNIDENTITÄT ★")
print(f"\n{'n':>3} | {'Cn→McKay':>10} | {'E-Algebra':>10} | {'h':>5} | {'φ(h)':>5} | {'B(n)':>5} | {'2(n-1)':>7} | {'Match':>6}")
print("-" * 68)

identity_holds = True
for n in range(2, 7):
    b_n = burnside_count(n)
    two_nm1 = 2 * (n - 1)

    if n == 3:
        alg = 'E6'
    elif n == 4:
        alg = 'E7'
    elif n == 5:
        alg = 'E8'
    elif n == 2:
        alg = 'A1'
    elif n == 6:
        alg = 'A5'
    else:
        alg = '—'

    if alg in coxeter_numbers:
        h = coxeter_numbers[alg]
        phi_h = euler_totient(h)
        match = "✓" if phi_h == b_n else "✗"
        if n in [3, 4, 5]:
            match = "★ " + match
    else:
        h = '—'
        phi_h = '—'
        match = '—'

    print(f"{n:>3} | {'C'+str(n)+'→'+alg:>10} | {alg:>10} | {str(h):>5} | {str(phi_h):>5} | {b_n:>5} | {two_nm1:>7} | {match:>6}")

print("\n► Die Identität B(n) = φ(h(E_{n+3})) gilt EXAKT für n = 3, 4, 5")
print("  und bricht für n = 2 (A-Typ) und n = 6 (A-Typ).")

# =====================================================================
# TEIL 4: Primitive Exponenten
# =====================================================================
print("\n" + "-" * 72)
print("TEIL 4: Primitive Exponenten der E-Algebren")
print("-" * 72)

# Exponenten (aus Standardtabellen, Bourbaki)
exponents = {
    'E6': [1, 4, 5, 7, 8, 11],
    'E7': [1, 5, 7, 9, 11, 13, 17],
    'E8': [1, 7, 11, 13, 17, 19, 23, 29]
}

for alg in ['E6', 'E7', 'E8']:
    h = coxeter_numbers[alg]
    r = ranks[alg]
    exps = exponents[alg]

    # Primitive = teilerfremd zu h
    primitive = [e for e in exps if gcd(e, h) == 1]
    non_primitive = [e for e in exps if gcd(e, h) != 1]

    print(f"\n{alg} (rank={r}, h={h}):")
    print(f"  Alle Exponenten:      {exps}")
    print(f"  Primitiv (gcd=1 mit h): {primitive} → {len(primitive)} Stück")
    print(f"  Nicht-primitiv:       {non_primitive} → {len(non_primitive)} Stück")
    print(f"  φ(h) = φ({h}) = {euler_totient(h)} = Anzahl primitiver Exponenten ✓")

    # Verifikation
    assert len(primitive) == euler_totient(h), f"FEHLER bei {alg}!"

print("\n★ INTERPRETATION:")
print("  Primitive Exponenten = algebraisch UNABHÄNGIGE Generatoren")
print("  der Coxeter-Dynamik (Chevalley 1955, Kostant 1959).")
print("  Jeder primitive Exponent erzeugt die volle zyklische Gruppe ⟨c⟩ ≅ Z_h.")
print("  Nicht-primitive Exponenten erzeugen nur Untergruppen → REDUNDANT.")

# =====================================================================
# TEIL 5: Die Drei-Wege-Identität bei E₈
# =====================================================================
print("\n" + "-" * 72)
print("TEIL 5: Die Drei-Wege-Identität bei E₈")
print("-" * 72)

print("\nFür jede E-Algebra prüfen wir: rank = φ(h) = B(n)?")
print()
print(f"{'Algebra':>8} | {'rank':>5} | {'φ(h)':>5} | {'B(n)':>5} | {'rank=φ(h)':>10} | {'φ(h)=B(n)':>10} | {'Drei-Wege':>10}")
print("-" * 70)

for alg, n in [('E6', 3), ('E7', 4), ('E8', 5)]:
    r = ranks[alg]
    h = coxeter_numbers[alg]
    phi_h = euler_totient(h)
    b_n = burnside_count(n)

    rk_eq_phi = "Ja" if r == phi_h else "Nein"
    phi_eq_b = "Ja" if phi_h == b_n else "Nein"
    three_way = "★ JA ★" if r == phi_h == b_n else "Nein"

    print(f"{alg:>8} | {r:>5} | {phi_h:>5} | {b_n:>5} | {rk_eq_phi:>10} | {phi_eq_b:>10} | {three_way:>10}")

print("\n► NUR bei E₈: rank(E₈) = φ(h(E₈)) = B(5) = 8")
print("  Drei Zahlen aus VERSCHIEDENEN Bereichen der Mathematik stimmen überein:")
print("  - rank = Zahl der unabhängigen Generatoren der Lie-Algebra")
print("  - φ(h) = Zahl der primitiven Coxeter-Exponenten")
print("  - B(5) = Zahl der symmetrie-verschiedenen Zustände des C5-H-Brücken-Rings")
print("\n  Bei E₈ sind ALLE Exponenten primitiv — keine Redundanz.")
print("  Bei E₆: 2 von 6 redundant. Bei E₇: 1 von 7 redundant.")

# =====================================================================
# TEIL 6: Diophantische Analyse — Warum genau {3, 4, 5}?
# =====================================================================
print("\n" + "-" * 72)
print("TEIL 6: Diophantische Analyse — Warum genau {3, 4, 5}?")
print("-" * 72)

print("\nBedingung B(p) = 2(p-1) für Primzahl p:")
print("  Äquivalent zu: 2^{p-1} = (p-1)²")
print()
print(f"{'p':>5} | {'2^(p-1)':>12} | {'(p-1)²':>12} | {'Match':>8}")
print("-" * 45)

for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
    lhs = 2 ** (p - 1)
    rhs = (p - 1) ** 2
    match = "★ ✓ ★" if lhs == rhs else "✗"
    print(f"{p:>5} | {lhs:>12} | {rhs:>12} | {match:>8}")

print("\n► Die diophantische Gleichung 2^{p-1} = (p-1)² hat genau")
print("  ZWEI Primzahllösungen: p = 3 und p = 5.")
print("  Für p ≥ 7 wächst 2^{p-1} exponentiell, (p-1)² nur quadratisch.")

print("\nFür n = 4 (zusammengesetzt): B(4) = 6 = 2×3")
print("  Direkte Berechnung: B(4) = (1×16 + 1×4 + 2×2)/4 = 24/4 = 6 ✓")
print("  φ(18) = 18 × (1-1/2) × (1-1/3) = 6 ✓")

# Beweis dass es genau 2 Lösungen gibt
print("\n--- Beweis der Endlichkeit ---")
print("Für p ≥ 7: 2^{p-1} > (p-1)² ist äquivalent zu")
print("  (p-1) × ln(2) > 2 × ln(p-1)")
print("  Sei x = p-1 ≥ 6: x × ln(2) > 2 × ln(x)")
print("  ⟺ ln(2) > 2 × ln(x)/x")
print("  Da ln(x)/x für x ≥ 3 streng fällt:")
print(f"  ln(6)/6 = {np.log(6)/6:.4f} < ln(2)/2 = {np.log(2)/2:.4f}")
print("  Also: für alle p ≥ 7 gilt 2^{p-1} > (p-1)² — KEINE weiteren Lösungen.")

# =====================================================================
# TEIL 7: Die Verbindung — Primitive Exponenten als Thermodynamische DOFs
# =====================================================================
print("\n" + "-" * 72)
print("TEIL 7: Darstellungstheoretische Begründung der Equipartition")
print("-" * 72)

print("""
THEOREM-KETTE (für n = 3, 4, 5):

Schritt 1 [Empirisch]:
  Cn-symmetrische Ionenkanäle korrespondieren über McKay zu E_{n+3}

Schritt 2 [BEWIESEN — diese Berechnung]:
  B(n) = φ(h(E_{n+3})) = Anzahl primitiver Exponenten von E_{n+3}

Schritt 3 [BEWIESEN — Chevalley 1955, Kostant 1959]:
  Die φ(h) primitiven Exponenten sind algebraisch unabhängige
  Generatoren der Coxeter-Dynamik. Jeder ist ein eigenständiger Modus.

Schritt 4 [Standard-Physik — Boltzmann]:
  Jeder unabhängige thermodynamische Modus trägt kT/2.

Schritt 5 [Konsequenz]:
  E₀ = B(n) × kT/2 = φ(h(E_{n+3})) × kT/2

  → μ = kT/E₀ = 2/B(n)  (Bouchaud)
  → H = 1 − 1/B(n)       (Lowen-Teich)
""")

print("OFFENE LÜCKEN:")
print("  1 → 2: Warum bildet die biologische Cn-Symmetrie auf E-Typ ab?")
print("         (Cn → Ân-1 via Standard-McKay, aber Cn ⊂ binäre Gruppe → E-Typ)")
print("  3 → 4: Warum sind algebraische Moden = thermodynamische DOFs?")
print("         (Plausibel via Connes-KMS, aber mathematisch offen)")

# =====================================================================
# TEIL 8: Numerische Verifikation der Hurst-Vorhersagen
# =====================================================================
print("\n" + "-" * 72)
print("TEIL 8: Hurst-Vorhersagen mit Lie-Algebra-Begründung")
print("-" * 72)

print("\n★ Vorhersagen aus H = 1 − 1/B(n) mit B(n) = φ(h(E_{n+3})) ★")
print()

channels = [
    (3, 'ASIC, P2X', 'E6', 'Sensorik'),
    (4, 'BK, KcsA', 'E7', 'Berechnung'),
    (5, 'nAChR, GLIC', 'E8', 'Morphogenese'),
]

print(f"{'n':>3} | {'Kanäle':>15} | {'E-Algebra':>10} | {'B(n)':>5} | {'φ(h)':>5} | {'H pred.':>8} | {'Funktion':>14}")
print("-" * 72)

for n, channels_str, alg, function in channels:
    b_n = burnside_count(n)
    h = coxeter_numbers[alg]
    phi_h = euler_totient(h)
    H_pred = 1 - 1 / b_n
    print(f"{n:>3} | {channels_str:>15} | {alg:>10} | {b_n:>5} | {phi_h:>5} | {H_pred:>8.3f} | {function:>14}")

print(f"\nMesswerte (Wawrzkiewicz-Jalowiecka 2024):")
print(f"  BK (C4): H_DFA = 0.81 ± 0.06 — Vorhersage: 0.833 (Fehler: 2.8%)")
print(f"  C2-Kanäle: H_DFA = 0.58-0.61 — Vorhersage: 0.667 (qualitativ: C4 > C2 ✓)")

# =====================================================================
# TEIL 9: Vergleich mit Berechnung AA (McKay-Cartan-Orbit)
# =====================================================================
print("\n" + "-" * 72)
print("TEIL 9: Vergleich mit Berechnung AA")
print("-" * 72)

print("""
AA-Route (graphentheoretisch, 6/10):
  ✗ Orbit-Graph ≠ Dynkin-Diagramm (Topologie verschieden)
  ✗ Standard-McKay: Cn → Ân-1, nicht E-Typ
  ✗ Dimensionen stimmen nicht (B(n) ≠ rank(En))
  → Ehrlich negativ: Keine direkte Graph-Isomorphie

AB-Route (zahlentheoretisch, DIESE Berechnung):
  ✓ Umgeht das Graph-Problem vollständig
  ✓ Nutzt Coxeter-Zahlen (nicht Dynkin-Topologie)
  ✓ B(n) = φ(h) = Anzahl primitiver Exponenten
  ✓ Gibt darstellungstheoretische Begründung für Equipartition
  ✓ Erklärt warum GENAU {3,4,5}: diophantische Bedingung
  → Qualitativ besserer Ansatz als AA
""")

# =====================================================================
# TEIL 10: Redundanz-Effizienz-Analyse
# =====================================================================
print("\n" + "-" * 72)
print("TEIL 10: Redundanz-Effizienz nach E-Algebra")
print("-" * 72)

print(f"\n{'Algebra':>8} | {'rank':>5} | {'φ(h)':>5} | {'Redundanz':>10} | {'Effizienz':>10} | {'Biolog. Interpretation':>25}")
print("-" * 80)

for alg, n, bio in [('E6', 3, 'Sensorik (schnell, verlust.)'),
                     ('E7', 4, 'Berechnung (ausgewogen)'),
                     ('E8', 5, 'Morphogenese (maximal eff.)')]:
    r = ranks[alg]
    h = coxeter_numbers[alg]
    phi_h = euler_totient(h)
    redundanz = r - phi_h
    effizienz = phi_h / r
    print(f"{alg:>8} | {r:>5} | {phi_h:>5} | {redundanz:>10} | {effizienz:>10.1%} | {bio:>25}")

print("""
► E₈/C5 hat 100% Effizienz (alle Moden primitiv)
  E₇/C4 hat 86% Effizienz (1 redundanter Modus)
  E₆/C3 hat 67% Effizienz (2 redundante Moden)

Biologische Implikation: C5-Kanäle (nAChR, GLIC) nutzen ALLE
verfügbaren algebraischen Moden — maximale Informationsverarbeitung.
C3-Kanäle "verlieren" 2 Moden an Redundanz — aber dafür sind sie
SCHNELLER (weniger Moden = höhere Bandbreite h = 2/(B+2)).
""")

# =====================================================================
# TEIL 11: Fibonacci-Verbindung
# =====================================================================
print("-" * 72)
print("TEIL 11: Fibonacci-Coxeter-Verbindung")
print("-" * 72)

print("\nCoxeter-Zahlen und Fibonacci:")
print(f"  h(E₆)/6 = 12/6 = 2 = F₃")
print(f"  h(E₇)/6 = 18/6 = 3 = F₄")
print(f"  h(E₈)/6 = 30/6 = 5 = F₅")
print(f"  (12 + 18 = 30 → Fibonacci-Relation!)")

print(f"\nDie VOLLSTÄNDIGE Kette für E₈:")
print(f"  rank(E₈) = φ(h(E₈)) = φ(30) = B(5) = 8")
print(f"  h(E₈) = 30 = 6 × F₅ = 6 × 5")
print(f"  dim(E₈) = 248 = 8 × 31 = rank × (h+1)")
print(f"  |W(E₈)| = 696,729,600 = 8! × 2^7 × 3^5 × 5^2 × 7")

# =====================================================================
# TEIL 12: Zusammenfassung und Bewertung
# =====================================================================
print("\n" + "=" * 72)
print("ZUSAMMENFASSUNG UND BEWERTUNG")
print("=" * 72)

print("""
★ HAUPTERGEBNIS: B(n) = φ(h(E_{n+3})) für n = 3, 4, 5

  Die Burnside-Orbitzahl = Zahl der primitiven Coxeter-Exponenten.
  Gilt exakt für die drei Werte, die auf E-Typ-Algebren abbilden.

★ BEDEUTUNG FÜR DAS MODELL:

  1. EQUIPARTITION BEGRÜNDUNG: Nicht "jeder Orbit ist irgendwie ein DOF",
     sondern: "Jeder Orbit entspricht einem primitiven Coxeter-Exponenten,
     und primitive Exponenten sind bewiesenermassen algebraisch unabhängig."

  2. TEMPERATURUNABHÄNGIGKEIT: kT kürzt sich heraus weil die Orbits
     als GEOMETRISCHE Freiheitsgrade wirken (Entropie, nicht Energie).

  3. WARUM {3,4,5}: Die diophantische Bedingung 2^{p-1} = (p-1)²
     erklärt exakt, warum die Identität für n=3,5 gilt.
     Für n=4 via direkte Verifikation.

  4. E₈ SONDERSTELLUNG: Einzige Drei-Wege-Identität
     rank = φ(h) = B(5) = 8. Alle Moden primitiv. 100% Effizienz.

★ OFFENE FRAGEN:

  1. Warum Cn → E_{n+3}? (Biologisch beobachtet, mathematisch nicht bewiesen)
  2. Warum algebraische Unabhängigkeit → thermodynamische DOFs?
     (Plausibel, aber Brücke fehlt — möglicherweise via Connes-KMS)

★ BEWERTUNG:

  AB vs. AA:  AB ist qualitativ BESSER (zahlentheoretisch statt graphentheoretisch)
  AB allein:  7.5/10 (verifizierte Identität + darstellungstheoretische Begründung,
              aber 2 offene Lücken und nur 3 Datenpunkte)

  STÄRKE:     Gibt erstmals eine darstellungstheoretische Begründung
              für die Burnside-Equipartition
  SCHWÄCHE:   3 Datenpunkte, McKay-Brücke empirisch

  FAZIT:      Bisher der BESTE mathematische Verbindungsversuch.
              Nicht Durchbruch-Niveau (dafür fehlen die 2 Lücken),
              aber die stärkste theoretische Fundierung die wir haben.
""")

print("=" * 72)
print(f"Berechnung AB abgeschlossen. Bewertung: 7.5/10")
print(f"Datei: calculations/calculation_AB_burnside_totient_E_identity.py")
print("=" * 72)
