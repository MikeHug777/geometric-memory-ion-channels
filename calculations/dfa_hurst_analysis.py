"""
DFA Hurst Analysis — Detrended Fluctuation Analysis for Ion Channel Gating
============================================================================
Author:  Michael Hug
Date:    2026-02-22
Project: Interface First — Geometric Memory
Paper:   Paper 2, Analysis Pipeline

Original contribution:
  DFA implementation optimized for binary single-channel time series,
  with scale selection (s = 10 to N/4) appropriate for ion channel
  recordings. Core analysis tool for extracting Hurst exponents from
  patch-clamp data as part of the Geometric Memory validation pipeline.

Dependencies: numpy, matplotlib, scipy
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def calc_rms(x, scale):
    """
    Calculate Root Mean Square (RMS) for a given scale (window size).
    """
    shape = (x.shape[0] // scale, scale)
    X = np.reshape(x[:shape[0] * shape[1]], shape)
    
    rms = np.zeros(shape[0])
    for i in range(shape[0]):
        # Fit a linear polynomial (degree 1) to detrend the window
        p = np.polyfit(np.arange(scale), X[i], 1)
        trend = np.polyval(p, np.arange(scale))
        # Calculate RMS of the detrended window
        rms[i] = np.sqrt(np.mean((X[i] - trend) ** 2))
    
    return np.mean(rms)

def dfa(time_series, scales=None):
    """
    Perform Detrended Fluctuation Analysis (DFA) on a time series.
    Returns scales, fluctuations (RMS), and the Hurst exponent (H).
    """
    # 1. Integrate the time series (cumulative sum after subtracting the mean)
    y = np.cumsum(time_series - np.mean(time_series))
    
    if scales is None:
        # Default scales: logarithmic spacing from 10 to N/4
        min_scale = 10
        max_scale = len(y) // 4
        scales = np.unique(np.logspace(np.log10(min_scale), np.log10(max_scale), 30).astype(int))
    
    # 2. Calculate RMS for each scale
    fluctuations = np.array([calc_rms(y, scale) for scale in scales])
    
    # 3. Fit a line in the log-log space to find the Hurst exponent (slope)
    # Exclude any NaNs or zeros that might have crept in
    valid = fluctuations > 0
    scales_valid = scales[valid]
    fluctuations_valid = fluctuations[valid]
    
    log_scales = np.log10(scales_valid)
    log_fluctuations = np.log10(fluctuations_valid)
    
    # Linear fit: log(F) = H * log(scale) + C
    p = np.polyfit(log_scales, log_fluctuations, 1)
    hurst_exponent = p[0]
    
    return scales_valid, fluctuations_valid, hurst_exponent

def analyze_and_plot(time_series, channel_name="Unknown Channel", expected_h=None):
    """
    Run DFA and plot the results.
    """
    scales, fluctuations, H = dfa(time_series)
    
    plt.figure(figsize=(8, 6))
    plt.loglog(scales, fluctuations, 'bo', label='DFA Fluctuation')
    
    # Plot the fit line
    fit_line = 10**(np.polyval(np.polyfit(np.log10(scales), np.log10(fluctuations), 1), np.log10(scales)))
    plt.loglog(scales, fit_line, 'r-', linewidth=2, label=f'Fit (H = {H:.3f})')
    
    if expected_h is not None:
        plt.title(f'DFA: {channel_name} (Calculated H={H:.3f}, Expected H={expected_h:.3f})')
    else:
        plt.title(f'DFA: {channel_name} (Calculated H={H:.3f})')
        
    plt.xlabel('Scale (log10)')
    plt.ylabel('Fluctuation RMS (log10)')
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.show()
    
    print(f"[{channel_name}] Calculated Hurst Exponent (H): {H:.3f}")
    if expected_h is not None:
        print(f"[{channel_name}] Expected Hurst Exponent (H): {expected_h:.3f}")
        error = abs(H - expected_h) / expected_h * 100
        print(f"[{channel_name}] Deviation from Theory: {error:.2f}%")
        
    return H

# Example Usage:
if __name__ == "__main__":
    # Simulate a time series (replace this with actual patch-clamp data: 0 for closed, 1 for open)
    # For demonstration, we use fractional Gaussian noise (requires fbm library, but here we just use random noise as a placeholder)
    # In reality, you'd load your binary open/closed sequence here.
    
    print("WARNING: This example uses white noise (H=0.5). Replace with real patch-clamp data.")
    np.random.seed(42)
    synthetic_data = np.random.randn(10000) 
    
    # Analyze the synthetic data
    analyze_and_plot(synthetic_data, channel_name="Synthetic White Noise", expected_h=0.500)
    
    # --- How to use with real data ---
    # 1. Load your patch clamp data (e.g., from a CSV or text file)
    # data = np.loadtxt('my_patch_clamp_data.csv')
    # 
    # 2. Ensure it's a binary sequence of states (e.g., 0=closed, 1=open) or dwell times.
    # 
    # 3. Run the analysis:
    # analyze_and_plot(data, channel_name="P2X7 (C3)", expected_h=0.750)
