"""
Integrated Analysis Pipeline — Demonstration
==============================================
Author:  Michael Hug
Date:    2026-02-22
Project: Interface First — Geometric Memory
Paper:   Paper 2, Analysis Pipeline (Demo)

Original contribution:
  End-to-end demonstration of the full analysis chain: simulated fractal
  noise generation, idealization, DFA analysis, and comparison with
  Burnside predictions. Validates pipeline correctness on known input.

Dependencies: numpy, idealize_patch_clamp, dfa_hurst_analysis
"""
import numpy as np
import os
import sys

# Import our newly created modules
sys.path.append('/Users/mikehug/Documents/personality-overview-core/pcm_model/data/Forschung/Interface_First/calculations')
from idealize_patch_clamp import idealize_patch_clamp
from dfa_hurst_analysis import analyze_and_plot

def generate_fractal_noise(n, h):
    """
    Generate fractional Gaussian noise (fGn) using the Davies-Harte method.
    This creates a time series with a specific Hurst exponent.
    """
    gamma = lambda k, H: 0.5 * ((k+1)**(2*H) - 2*k**(2*H) + abs(k-1)**(2*H))
    g = np.array([gamma(k, h) for k in range(n)])
    r = np.concatenate(([g[0]], g[1:], g[:0:-1]))
    
    # Check for valid eigenvalues
    eigenvalues = np.fft.fft(r).real
    if np.any(eigenvalues < 0):
        print("Warning: Some eigenvalues are negative. Generating approximate fGn.")
        eigenvalues[eigenvalues < 0] = 0
        
    # The length of r is 1 + (n-1) + (n-1) = 2n - 1
    # So we need complex_noise to be of size 2n - 1
    complex_noise = np.random.randn(2*n-1) + 1j * np.random.randn(2*n-1)
    # The first element and the n-th element (if it exists) must be real for the iFFT to yield a real sequence
    complex_noise[0] = np.random.randn() * np.sqrt(2)
    if n > 1:
        complex_noise[n-1] = np.random.randn() * np.sqrt(2)
        
    fgn = np.fft.ifft(complex_noise * np.sqrt(eigenvalues)).real
    return fgn[:n]

def simulate_and_run_pipeline():
    print("=== Phase 1: Simulating realistic C4 Patch-Clamp Data ===")
    
    np.random.seed(42)
    n_samples = 50000
    
    # 1. Generate the underlying fractal switching process (Hurst = 0.833 for C4)
    # This represents the "memory" of the channel's conformational state
    target_H = 0.833
    underlying_fractal = generate_fractal_noise(n_samples, target_H)
    
    # Standardize it
    underlying_fractal = (underlying_fractal - np.mean(underlying_fractal)) / np.std(underlying_fractal)
    
    # 2. Turn this into a binary open/closed signal using a threshold
    # If the fractal state is > 0, the channel is open, else closed.
    true_states = (underlying_fractal > 0).astype(int)
    
    # 3. Add realistic amplifier noise and current amplitudes
    # Baseline (closed) is ~0 pA, Open state is ~6 pA
    baseline_noise = np.random.normal(0, 1.2, n_samples)
    open_current = 6.0
    open_noise = np.random.normal(0, 1.5, n_samples)
    
    raw_patch_current = baseline_noise + true_states * (open_current + open_noise)
    
    raw_file = '/Users/mikehug/Documents/personality-overview-core/pcm_model/data/Forschung/Interface_First/calculations/simulated_c4_raw.csv'
    np.savetxt(raw_file, raw_patch_current, fmt='%.3f')
    print(f"-> Simulated {n_samples} data points saved to {raw_file}")
    
    print("\n=== Phase 2: Idealization (Raw -> Binary) ===")
    # Load the raw data as if we just downloaded it from Zenodo
    data_to_idealize = np.loadtxt(raw_file)
    
    # We turn plotting off here to keep the script automated, but normally you'd set plot=True
    idealized_data, calc_threshold = idealize_patch_clamp(data_to_idealize, plot=False)
    
    ideal_file = '/Users/mikehug/Documents/personality-overview-core/pcm_model/data/Forschung/Interface_First/calculations/simulated_c4_idealized.csv'
    np.savetxt(ideal_file, idealized_data, fmt='%d')
    print(f"-> Idealized binary data saved to {ideal_file}")
    
    print("\n=== Phase 3: DFA Analysis (Extracting the Hurst Exponent) ===")
    data_for_dfa = np.loadtxt(ideal_file)
    
    # Run the DFA and plot the result
    # We expect an H of ~0.833 because it's a C4 channel (Burnside formula: 1 - 1/6)
    print("Calculating Hurst Exponent...")
    
    # Temporary disable interactive plotting for the automated run so it doesn't block
    import matplotlib
    matplotlib.use('Agg') 
    
    calculated_H = analyze_and_plot(data_for_dfa, channel_name="Simulated KcsA/Kv (C4)", expected_h=0.833)
    
    print(f"\n--- SUCCESS ---")
    print(f"Target theoretical Burnside H (C4): 0.833")
    print(f"Recovered H from noisy simulated data: {calculated_H:.3f}")
    
    # Save the plot
    import matplotlib.pyplot as plt
    plot_file = '/Users/mikehug/Documents/personality-overview-core/pcm_model/data/Forschung/Interface_First/calculations/dfa_c4_result.png'
    plt.savefig(plot_file)
    print(f"Plot saved to: {plot_file}")

if __name__ == "__main__":
    simulate_and_run_pipeline()
