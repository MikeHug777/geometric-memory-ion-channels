#!/usr/bin/env python3
"""
Calculation AH: Relaxation Spectrum of Cn-Symmetric Trap Models
================================================================
Author:  Michael Hug
Date:    2026-02-25
Project: Interface First — Geometric Memory
Paper:   Paper 3, SI Section

Key question: Do the Coxeter invariant degrees of E_{n+3} appear
as relaxation timescales in a Bouchaud trap model on the Cn-symmetric
binary ring {0,1}^n?

The Cn rotation groups the 2^n binary states into B(n) Burnside orbits.
The hypercube graph (single-bit-flip transitions) induces an orbit
adjacency structure. We add Arrhenius trap depths and compute the
eigenvalue spectrum of the resulting rate matrix.

Tests:
  TEST 1: Orbit structure — enumerate orbits, adjacency matrix, sizes
  TEST 2: Topological eigenvalues — bare spectrum without trap depths
  TEST 3: Arrhenius eigenvalues — spectrum with Exp(E₀) trap depths
  TEST 4: Degree-ratio comparison — τ_j/τ_1 vs d_j/d_1 for E-series
  TEST 5: MC autocorrelation — verify power-law C(τ) ~ τ^{-1/B}

Dependencies: numpy, matplotlib, scipy
"""

import numpy as np
import matplotlib.pyplot as plt
from math import gcd
import os
import warnings
warnings.filterwarnings('ignore')

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "mc_relaxation_spectrum")
os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# COXETER DATA
# ============================================================

E_SERIES = {
    3: {  # C3 → E6
        'algebra': 'E6', 'rank': 6, 'h': 12,
        'degrees': [2, 5, 6, 8, 9, 12],
        'exponents': [1, 4, 5, 7, 8, 11],
        'prim_degrees': [2, 6, 8, 12],  # gcd(e, 12) = 1
        'color': '#2196F3',
    },
    4: {  # C4 → E7
        'algebra': 'E7', 'rank': 7, 'h': 18,
        'degrees': [2, 6, 8, 10, 12, 14, 18],
        'exponents': [1, 5, 7, 9, 11, 13, 17],
        'prim_degrees': [2, 6, 8, 12, 14, 18],  # gcd(e, 18) = 1
        'color': '#4CAF50',
    },
    5: {  # C5 → E8
        'algebra': 'E8', 'rank': 8, 'h': 30,
        'degrees': [2, 8, 12, 14, 18, 20, 24, 30],
        'exponents': [1, 7, 11, 13, 17, 19, 23, 29],
        'prim_degrees': [2, 8, 12, 14, 18, 20, 24, 30],  # all primitive
        'color': '#F44336',
    },
}


# ============================================================
# HELPERS
# ============================================================

def euler_phi(m):
    result = m
    p = 2
    temp = m
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result


def burnside_number(n):
    total = 0
    for d in range(1, n + 1):
        if n % d == 0:
            total += euler_phi(d) * (2 ** (n // d))
    return total // n


def cyclic_rotate(state, n):
    """Rotate binary string of length n by one position."""
    return ((state >> 1) | ((state & 1) << (n - 1))) & ((1 << n) - 1)


# ============================================================
# TEST 1: Orbit Structure
# ============================================================

def enumerate_orbits(n):
    """Enumerate all Cn orbits on {0,1}^n.

    Returns:
        orbits: list of lists, each inner list contains the states in one orbit
        state_to_orbit: dict mapping state -> orbit index
    """
    seen = set()
    orbits = []
    state_to_orbit = {}

    for state in range(2 ** n):
        if state in seen:
            continue
        orbit = []
        s = state
        for _ in range(n):
            if s not in seen:
                orbit.append(s)
                seen.add(s)
            s = cyclic_rotate(s, n)
        orbits.append(sorted(orbit))

    # Build state-to-orbit map
    for i, orbit in enumerate(orbits):
        for state in orbit:
            state_to_orbit[state] = i

    return orbits, state_to_orbit


def orbit_adjacency(orbits, state_to_orbit, n):
    """Compute the orbit adjacency matrix from the hypercube graph.

    A[i,j] = average number of neighbors in orbit j for a state in orbit i.
    Neighbors are states differing by exactly one bit (single-bit flip).
    """
    B = len(orbits)
    A = np.zeros((B, B))
    orbit_sizes = np.array([len(orbit) for orbit in orbits])

    for i, orbit in enumerate(orbits):
        for state in orbit:
            for bit in range(n):
                neighbor = state ^ (1 << bit)
                j = state_to_orbit[neighbor]
                A[i, j] += 1
        # Average over states in orbit i
        A[i] /= orbit_sizes[i]

    return A, orbit_sizes


def test1_orbit_structure():
    """Enumerate orbits and compute adjacency for n = 3, 4, 5, 6."""
    print("=" * 72)
    print("TEST 1: Orbit Structure of Cn on {0,1}^n")
    print("=" * 72)

    results = {}

    for n in [3, 4, 5, 6]:
        orbits, s2o = enumerate_orbits(n)
        A, sizes = orbit_adjacency(orbits, s2o, n)
        B = len(orbits)

        print(f"\n  C{n}: 2^{n} = {2**n} states → B({n}) = {B} orbits")
        print(f"  {'Orbit':<6} {'Size':<6} {'States':<40} {'Self-loops':<10} {'Neighbors'}")
        print("  " + "-" * 75)

        for i in range(B):
            states_str = str(orbits[i][:6])
            if len(orbits[i]) > 6:
                states_str = states_str[:-1] + ", ...]"
            self_loops = A[i, i]
            neighbors = sum(A[i, j] for j in range(B) if j != i)
            print(f"  {i:<6} {sizes[i]:<6} {states_str:<40} "
                  f"{self_loops:<10.1f} {neighbors:.1f}")

        results[n] = {
            'orbits': orbits, 'state_to_orbit': s2o,
            'adjacency': A, 'sizes': sizes, 'B': B
        }

    return results


# ============================================================
# TEST 2: Topological Eigenvalues (no trap depths)
# ============================================================

def test2_topological_eigenvalues(orbit_data):
    """Eigenvalues of the orbit transition matrix (bare, no trap depths).

    The transition matrix P[i,j] = A[i,j] / n gives the probability of
    jumping from orbit i to orbit j via a random single-bit flip.
    """
    print("\n" + "=" * 72)
    print("TEST 2: Topological Eigenvalues (no trap depths)")
    print("=" * 72)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    results = {}

    for idx, n in enumerate([3, 4, 5, 6]):
        A = orbit_data[n]['adjacency']
        B = orbit_data[n]['B']

        # Transition matrix: P = A / n (each state has n neighbors)
        P = A / n

        # Eigenvalues
        eigenvalues = np.sort(np.linalg.eigvals(P).real)[::-1]

        # Relaxation times: τ_j = -1/ln(|λ_j|) for discrete time
        # or τ_j = 1/(1 - λ_j) for the continuous-time version
        relax_times = []
        for lam in eigenvalues:
            if abs(1 - lam) > 1e-10:  # skip equilibrium eigenvalue
                relax_times.append(1.0 / (1.0 - lam))

        relax_times = sorted(relax_times, reverse=True)

        print(f"\n  C{n} (B={B}):")
        print(f"    Eigenvalues of P: {np.array2string(eigenvalues, precision=4)}")
        print(f"    Relaxation times τ = 1/(1-λ):")
        for j, tau in enumerate(relax_times):
            ratio = tau / relax_times[-1] if relax_times[-1] > 0 else 0
            print(f"      τ_{j+1} = {tau:.4f}  (ratio to fastest: {ratio:.3f})")

        # Compare with Coxeter degrees if available
        if n in E_SERIES:
            prim_deg = E_SERIES[n]['prim_degrees']
            print(f"    Coxeter primitive degree ratios (d/d_min):")
            d_min = min(prim_deg)
            for d in prim_deg:
                print(f"      d = {d}, d/d_min = {d/d_min:.1f}")

        results[n] = {
            'eigenvalues': eigenvalues,
            'relax_times': relax_times,
        }

        # Plot
        ax = axes[idx // 2][idx % 2]
        ax.stem(range(B), eigenvalues, linefmt='-', markerfmt='o',
                basefmt='k-')
        ax.axhline(0, color='gray', ls=':', alpha=0.5)
        ax.axhline(1, color='gray', ls=':', alpha=0.5)

        # Mark hypercube eigenvalues
        hyp_eigs = [1 - 2 * k / n for k in range(n + 1)]
        for he in set(hyp_eigs):
            ax.axhline(he, color='red', ls='--', alpha=0.2, lw=0.5)

        ax.set_xlabel('Eigenvalue index')
        ax.set_ylabel('λ')
        ax.set_title(f'C{n}: B={B} orbits, {B} eigenvalues')
        ax.grid(True, alpha=0.3)

    plt.suptitle('Topological eigenvalues of the orbit transition matrix\n'
                 '(no trap depths — pure combinatorial structure)',
                 fontsize=12)
    plt.tight_layout()
    path = os.path.join(OUTDIR, "fig1_topological_eigenvalues.png")
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.savefig(path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"\n  → Saved: {path}")
    plt.close()

    return results


# ============================================================
# TEST 3: Arrhenius Eigenvalues (with trap depths)
# ============================================================

def test3_arrhenius_eigenvalues(orbit_data, n_realizations=500):
    """Eigenvalue spectrum of rate matrix with Exp(E₀) trap depths.

    Rate matrix: R[i,j] = A[i,j]/n * exp(-E_i/kT) for i ≠ j
                 R[i,i] = -Σ_{j≠i} R[i,j]

    The eigenvalues give relaxation rates. We average the LOG of the
    eigenvalue ratios over many realizations.
    """
    print("\n" + "=" * 72)
    print(f"TEST 3: Arrhenius Eigenvalues ({n_realizations} realizations)")
    print("=" * 72)

    kT = 1.0
    results = {}

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    for idx, n in enumerate([3, 4, 5, 6]):
        A = orbit_data[n]['adjacency']
        B = orbit_data[n]['B']
        E0 = burnside_number(n) * kT / 2

        all_log_ratios = []  # log(τ_j / τ_min) for each realization

        rng = np.random.default_rng(seed=42 + n)

        for r in range(n_realizations):
            # Draw trap depths
            trap_depths = rng.exponential(scale=E0, size=B)

            # Construct rate matrix
            R = np.zeros((B, B))
            for i in range(B):
                rate_i = np.exp(-trap_depths[i] / kT)
                for j in range(B):
                    if i != j:
                        R[i, j] = (A[i, j] / n) * rate_i
                R[i, i] = -np.sum(R[i, :])

            # Eigenvalues
            eigs = np.sort(np.linalg.eigvals(R).real)  # ascending
            # eigs[B-1] ≈ 0 (equilibrium), eigs[0] is most negative

            # Non-zero eigenvalues (skip the one closest to 0)
            nonzero = eigs[:-1]  # all except the largest (≈0)
            rates = -nonzero  # relaxation rates (positive)
            rates = rates[rates > 1e-30]  # filter numerical zeros

            if len(rates) >= 2:
                rates = np.sort(rates)  # ascending: slowest first
                tau = 1.0 / rates  # relaxation times: longest first
                log_ratios = np.log10(tau / tau[-1])  # ratio to fastest
                all_log_ratios.append(log_ratios)

        # Statistics over realizations
        min_len = min(len(lr) for lr in all_log_ratios)
        log_ratio_arr = np.array([lr[:min_len] for lr in all_log_ratios])
        median_log_ratios = np.median(log_ratio_arr, axis=0)
        q25 = np.percentile(log_ratio_arr, 25, axis=0)
        q75 = np.percentile(log_ratio_arr, 75, axis=0)

        print(f"\n  C{n} (B={B}, E₀={E0:.1f}kT):")
        print(f"    {min_len} relaxation times per realization")
        print(f"    Median log₁₀(τ_j/τ_min):")
        for j in range(min_len):
            print(f"      Mode {j+1}: {median_log_ratios[j]:.2f} "
                  f"[IQR: {q25[j]:.2f} – {q75[j]:.2f}]")

        # Compare with Coxeter degree ratios
        if n in E_SERIES:
            prim_deg = sorted(E_SERIES[n]['prim_degrees'])
            d_min = prim_deg[0]
            coxeter_log_ratios = np.log10(np.array(prim_deg) / d_min)
            print(f"    Coxeter log₁₀(d_j/d_min) for comparison:")
            for j, d in enumerate(prim_deg):
                print(f"      d={d}: {coxeter_log_ratios[j]:.3f}")

        results[n] = {
            'median_log_ratios': median_log_ratios,
            'q25': q25, 'q75': q75,
            'min_len': min_len,
            'log_ratio_arr': log_ratio_arr,
        }

        # Plot: distribution of log-relaxation-time ratios
        ax = axes[idx // 2][idx % 2]

        # Box plot of each mode's log-ratio
        positions = range(1, min_len + 1)
        bp = ax.boxplot([log_ratio_arr[:, j] for j in range(min_len)],
                        positions=positions, widths=0.6,
                        patch_artist=True,
                        medianprops=dict(color='black', lw=2))
        for patch in bp['boxes']:
            patch.set_facecolor(E_SERIES.get(n, {'color': '#9C27B0'})['color'])
            patch.set_alpha(0.4)

        # Overlay Coxeter degree ratios
        if n in E_SERIES:
            prim_deg = sorted(E_SERIES[n]['prim_degrees'])
            d_min = prim_deg[0]
            coxeter_lr = np.log10(np.array(prim_deg) / d_min)
            # Only plot if number matches
            if len(coxeter_lr) <= min_len + 1:
                for j, clr in enumerate(coxeter_lr):
                    if j < min_len:
                        ax.plot(j + 1, clr, 'r*', markersize=15, zorder=10,
                                label='Coxeter d_j/d₁' if j == 0 else '')

        ax.set_xlabel('Mode index (slowest → fastest)')
        ax.set_ylabel('log₁₀(τ_j / τ_min)')
        ax.set_title(f'C{n}: B={B}, {n_realizations} realizations')
        ax.grid(True, alpha=0.3)
        if n in E_SERIES:
            ax.legend(fontsize=9)

    plt.suptitle('Relaxation time ratios: MC (box) vs Coxeter degrees (★)\n'
                 'Rate matrix eigenvalues with Exp(E₀) trap depths',
                 fontsize=12)
    plt.tight_layout()
    path = os.path.join(OUTDIR, "fig2_arrhenius_eigenvalues.png")
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.savefig(path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"\n  → Saved: {path}")
    plt.close()

    return results


# ============================================================
# TEST 4: Degree-Ratio Comparison (scatter plot)
# ============================================================

def test4_degree_comparison(arrhenius_results):
    """Direct comparison: median τ_j/τ_min vs d_j/d_min."""
    print("\n" + "=" * 72)
    print("TEST 4: Relaxation Time Ratios vs Coxeter Degree Ratios")
    print("=" * 72)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for idx, n in enumerate([3, 4, 5]):
        ax = axes[idx]
        data = E_SERIES[n]
        ar = arrhenius_results[n]

        prim_deg = sorted(data['prim_degrees'])
        d_min = prim_deg[0]
        coxeter_ratios = np.array(prim_deg) / d_min

        # Median relaxation time ratios
        median_lr = ar['median_log_ratios']
        mc_ratios = 10 ** median_lr  # convert from log

        # Align: both should have B(n) - 1 or B(n) entries
        # Coxeter has φ(h) = B(n) primitive degrees
        # MC has B(n) - 1 relaxation times (one eigenvalue is 0)
        n_compare = min(len(coxeter_ratios), len(mc_ratios))

        if n_compare >= 2:
            # Sort both in ascending order
            cr = sorted(coxeter_ratios[:n_compare])
            mr = sorted(mc_ratios[:n_compare])

            ax.scatter(cr, mr, s=80,
                       c=data['color'], edgecolors='black', zorder=5)

            # Perfect match line
            max_val = max(max(cr), max(mr)) * 1.1
            ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3,
                    label='Perfect match')

            # Correlation
            if len(cr) >= 3:
                corr = np.corrcoef(cr, mr)[0, 1]
                ax.text(0.05, 0.95, f'r = {corr:.3f}',
                        transform=ax.transAxes, fontsize=12,
                        verticalalignment='top')

            print(f"\n  C{n} ({data['algebra']}):")
            print(f"    {'Coxeter d/d₁':<15} {'MC τ/τ_min':<15} {'Match?'}")
            print("    " + "-" * 40)
            for c, m in zip(cr, mr):
                match = "≈" if abs(c - m) / max(c, 0.01) < 0.3 else "≠"
                print(f"    {c:<15.3f} {m:<15.3f} {match}")

        ax.set_xlabel('Coxeter degree ratio d_j/d₁')
        ax.set_ylabel('MC relaxation time ratio τ_j/τ_min')
        ax.set_title(f'{data["algebra"]} (C{n}, B={burnside_number(n)})')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')

    plt.suptitle('Do Coxeter degrees predict relaxation timescales?',
                 fontsize=13)
    plt.tight_layout()
    path = os.path.join(OUTDIR, "fig3_degree_comparison.png")
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.savefig(path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"\n  → Saved: {path}")
    plt.close()


# ============================================================
# TEST 5: MC Autocorrelation
# ============================================================

def test5_mc_autocorrelation(orbit_data):
    """Simulate CTRW on orbit graph and compute autocorrelation.

    The autocorrelation C(τ) of the state signal should decay as
    a power law C(τ) ~ τ^{-(1-H)} = τ^{-1/B} for the Bouchaud model.

    We also compute the power spectrum P(ω) to look for discrete peaks.
    """
    print("\n" + "=" * 72)
    print("TEST 5: MC Autocorrelation and Power Spectrum")
    print("=" * 72)

    kT = 1.0

    # Parameters per symmetry
    mc_params = {
        3: {'T_max': 1e5,  'dt': 0.1, 'n_runs': 20},
        4: {'T_max': 1e6,  'dt': 1.0, 'n_runs': 10},
        5: {'T_max': 1e7,  'dt': 10.0, 'n_runs': 5},
    }

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    for idx, n in enumerate([3, 4, 5]):
        B_n = burnside_number(n)
        E0 = B_n * kT / 2
        mu = 2.0 / B_n
        params = mc_params[n]

        orbits = orbit_data[n]['orbits']
        A = orbit_data[n]['adjacency']
        B = orbit_data[n]['B']

        # Transition probabilities (normalized rows, excluding self-loops)
        P_trans = A.copy()
        for i in range(B):
            P_trans[i, i] = 0
            row_sum = P_trans[i].sum()
            if row_sum > 0:
                P_trans[i] /= row_sum

        print(f"\n  C{n} (B={B}, μ={mu:.4f}, {params['n_runs']} runs):")

        all_autocorr = []
        all_psd = []

        for run in range(params['n_runs']):
            rng = np.random.default_rng(seed=7000 + n * 1000 + run)

            # Draw trap depths for this realization
            trap_depths = rng.exponential(scale=E0, size=B)
            mean_wait = np.exp(trap_depths / kT)

            # Simulate CTRW
            T_max = params['T_max']
            dt = params['dt']
            n_samples = int(T_max / dt)

            # Event-based simulation
            current = rng.integers(B)
            t = 0.0
            events_t = [0.0]
            events_s = [current]

            while t < T_max:
                # Wait: exponential with mean = exp(E_current/kT)
                wait = rng.exponential(mean_wait[current])
                t += wait
                if t >= T_max:
                    break

                # Jump according to transition probs
                probs = P_trans[current]
                if probs.sum() > 0:
                    current = rng.choice(B, p=probs)
                else:
                    current = rng.integers(B)

                events_t.append(t)
                events_s.append(current)

            events_t = np.array(events_t)
            events_s = np.array(events_s)

            n_events = len(events_t) - 1
            if n_events < 100:
                print(f"    Run {run}: only {n_events} events — skipping")
                continue

            # Sample at regular intervals
            t_sample = np.arange(n_samples) * dt
            state_idx = np.searchsorted(events_t, t_sample, side='right') - 1
            state_idx = np.clip(state_idx, 0, len(events_s) - 1)
            signal = events_s[state_idx].astype(float)
            signal -= signal.mean()

            # Autocorrelation via FFT
            N = len(signal)
            F = np.fft.rfft(signal, n=2 * N)
            S = np.abs(F) ** 2
            C = np.fft.irfft(S)[:N] / N
            if C[0] > 0:
                C /= C[0]
            all_autocorr.append(C)

            # Power spectral density
            freq = np.fft.rfftfreq(N, d=dt)
            psd = np.abs(F[:len(freq)]) ** 2 * dt / N
            all_psd.append(psd)

        if not all_autocorr:
            print(f"    No valid runs for C{n}")
            continue

        # Average autocorrelation
        min_len_c = min(len(c) for c in all_autocorr)
        avg_C = np.mean([c[:min_len_c] for c in all_autocorr], axis=0)
        lags = np.arange(min_len_c) * params['dt']

        # Average PSD
        min_len_p = min(len(p) for p in all_psd)
        avg_psd = np.mean([p[:min_len_p] for p in all_psd], axis=0)
        freqs = freq[:min_len_p]

        # --- Plot autocorrelation ---
        ax = axes[0, idx]
        # Subsample for plotting
        n_plot = min(min_len_c, 10000)
        step = max(1, min_len_c // n_plot)
        plot_lags = lags[1:n_plot * step:step]
        plot_C = avg_C[1:n_plot * step:step]

        valid = plot_C > 0
        if valid.any():
            ax.loglog(plot_lags[valid], plot_C[valid], '.',
                      color=E_SERIES[n]['color'], markersize=2, alpha=0.5)

            # Fit power law in middle range
            fit_lo = len(plot_lags[valid]) // 5
            fit_hi = 4 * len(plot_lags[valid]) // 5
            if fit_hi > fit_lo + 10:
                x = np.log10(plot_lags[valid][fit_lo:fit_hi])
                y = np.log10(plot_C[valid][fit_lo:fit_hi])
                slope, intercept = np.polyfit(x, y, 1)
                exponent = -slope

                t_fit = np.logspace(np.log10(plot_lags[valid][fit_lo]),
                                    np.log10(plot_lags[valid][fit_hi]), 50)
                ax.loglog(t_fit, 10 ** intercept * t_fit ** slope,
                          'k--', lw=2,
                          label=f'slope = {slope:.3f}\n'
                                f'(theory: {-1/B_n:.3f})')
                ax.legend(fontsize=9)

                print(f"    C(τ) power-law exponent: {exponent:.4f} "
                      f"(theory: 1/B = {1/B_n:.4f})")

        ax.set_xlabel('Lag τ')
        ax.set_ylabel('C(τ)')
        ax.set_title(f'C{n}: Autocorrelation (B={B_n})')
        ax.grid(True, alpha=0.3, which='both')

        # --- Plot PSD ---
        ax2 = axes[1, idx]
        # Subsample
        valid_f = (freqs > 0) & (avg_psd > 0)
        if valid_f.any():
            step_f = max(1, valid_f.sum() // 3000)
            ax2.loglog(freqs[valid_f][::step_f], avg_psd[valid_f][::step_f],
                       '.', color=E_SERIES[n]['color'], markersize=1,
                       alpha=0.3)

            # Fit 1/f^α
            f_valid = freqs[valid_f]
            p_valid = avg_psd[valid_f]
            fit_lo_f = len(f_valid) // 10
            fit_hi_f = len(f_valid) // 2
            if fit_hi_f > fit_lo_f + 10:
                xf = np.log10(f_valid[fit_lo_f:fit_hi_f])
                yf = np.log10(p_valid[fit_lo_f:fit_hi_f])
                slope_f, intercept_f = np.polyfit(xf, yf, 1)

                f_fit = np.logspace(np.log10(f_valid[fit_lo_f]),
                                    np.log10(f_valid[fit_hi_f]), 50)
                ax2.loglog(f_fit,
                           10 ** intercept_f * f_fit ** slope_f,
                           'k--', lw=2,
                           label=f'slope = {slope_f:.2f}\n'
                                 f'(1/f^{{1-μ}}: {-(1-mu):.2f})')
                ax2.legend(fontsize=9)

                print(f"    P(ω) slope: {slope_f:.4f} "
                      f"(theory: {-(1-mu):.4f})")

            # Mark Coxeter degree frequencies (if meaningful)
            if n in E_SERIES:
                prim_deg = E_SERIES[n]['prim_degrees']
                for d in prim_deg:
                    f_coxeter = 1.0 / (d * params['dt'] * 10)
                    if f_valid[0] < f_coxeter < f_valid[-1]:
                        ax2.axvline(f_coxeter, color='red', alpha=0.15,
                                    lw=0.5)

        ax2.set_xlabel('Frequency ω')
        ax2.set_ylabel('P(ω)')
        ax2.set_title(f'C{n}: Power Spectrum')
        ax2.grid(True, alpha=0.3, which='both')

    plt.suptitle('MC: Autocorrelation C(τ) and Power Spectrum P(ω)\n'
                 'Bouchaud CTRW on Cn orbit graph',
                 fontsize=12)
    plt.tight_layout()
    path = os.path.join(OUTDIR, "fig4_mc_autocorrelation.png")
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.savefig(path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"\n  → Saved: {path}")
    plt.close()


# ============================================================
# SUMMARY
# ============================================================

def summary(topo_results, arrhenius_results):
    """Print final summary comparing all results."""
    print("\n" + "=" * 72)
    print("SUMMARY: Do Coxeter Degrees Appear as Relaxation Timescales?")
    print("=" * 72)

    print("""
  TEST 1 (Orbit Structure):
    Cn-orbits on {0,1}^n have well-defined adjacency from the hypercube.
    The orbit graph is non-trivial (not fully connected, not regular).

  TEST 2 (Topological Eigenvalues):
    The bare orbit transition matrix has eigenvalues drawn from the
    hypercube spectrum {1 - 2k/n : k = 0,...,n}.
    These do NOT match Coxeter degree ratios — the bare combinatorial
    structure of the binary ring does not encode E-algebra structure.

  TEST 3 (Arrhenius Eigenvalues):
    With Exp(E₀) trap depths, the rate matrix eigenvalues span many
    orders of magnitude. The ratios τ_j/τ_min vary strongly between
    realizations (heavy tails → wide distribution of timescales).
    """)

    for n in [3, 4, 5]:
        B = burnside_number(n)
        data = E_SERIES[n]
        ar = arrhenius_results[n]

        prim_deg = sorted(data['prim_degrees'])
        d_min = prim_deg[0]
        coxeter_ratios = np.log10(np.array(prim_deg) / d_min)

        median_lr = ar['median_log_ratios']
        n_compare = min(len(coxeter_ratios), len(median_lr))

        if n_compare >= 3:
            cr = coxeter_ratios[:n_compare]
            mr = median_lr[:n_compare]
            corr = np.corrcoef(cr, mr)[0, 1] if n_compare >= 3 else float('nan')
        else:
            corr = float('nan')

        print(f"  {data['algebra']} (C{n}, B={B}): "
              f"Correlation(log τ-ratios, log d-ratios) = {corr:.3f}")

    print("""
  INTERPRETATION:
    The NUMBER of effective modes (B(n) = φ(h)) determines the
    power-law exponent μ = 2/B and hence H = 1 - 1/B.
    Whether the SPECIFIC mode frequencies match Coxeter degrees
    is a separate, deeper question.

    Positive correlation → Coxeter degrees have physical meaning
                           beyond the mode count.
    No correlation      → Only B(n) = φ(h) matters, not the
                           individual degrees. The Molien partition
                           function captures the count correctly
                           but not the spectral structure.
    """)


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print()
    print("=" * 72)
    print("  Calculation AH: Relaxation Spectrum of Cn-Symmetric Trap Models")
    print("=" * 72)
    print()

    # Reference table
    print("Reference:")
    print(f"  {'n':>3} | {'B(n)':>5} | {'φ(h)':>5} | {'E₀/kT':>6} | "
          f"{'μ=2/B':>6} | {'H=1-1/B':>8} | Algebra")
    print("  " + "-" * 60)
    for n in [3, 4, 5, 6]:
        B = burnside_number(n)
        alg = E_SERIES[n]['algebra'] if n in E_SERIES else '(none)'
        h = E_SERIES[n]['h'] if n in E_SERIES else '—'
        phi = euler_phi(h) if isinstance(h, int) else '—'
        print(f"  {n:3d} | {B:5d} | {str(phi):>5} | {B/2:6.1f} | "
              f"{2/B:6.4f} | {1-1/B:8.4f} | {alg}")

    # Tests
    orbit_data = test1_orbit_structure()
    topo_results = test2_topological_eigenvalues(orbit_data)
    arrhenius_results = test3_arrhenius_eigenvalues(orbit_data)
    test4_degree_comparison(arrhenius_results)
    test5_mc_autocorrelation(orbit_data)

    # Summary
    summary(topo_results, arrhenius_results)

    print(f"\nAll figures saved to: {OUTDIR}/")
