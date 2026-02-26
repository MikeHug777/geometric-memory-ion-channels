"""
Idealize Patch Clamp — Binary State Detection from Raw Current Traces
======================================================================
Author:  Michael Hug
Date:    2026-02-22
Project: Interface First — Geometric Memory
Paper:   Paper 2, Analysis Pipeline

Original contribution:
  Half-amplitude threshold idealization converting raw patch-clamp
  current traces to binary open/closed time series. Preprocessing
  step for DFA analysis of single-channel fractal gating.

Dependencies: numpy, matplotlib, scipy
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def idealize_patch_clamp(raw_data, threshold=None, plot=True):
    """
    Idealize raw patch-clamp current data into binary open (1) and closed (0) states.
    
    Parameters:
    - raw_data: 1D numpy array of current values (e.g., in pA).
    - threshold: Cutoff value. If None, it will be automatically estimated using a histogram peak-finding approach.
    - plot: If True, plots a segment of the raw data vs. the idealized states.
    
    Returns:
    - idealized_states: 1D numpy array of 0s and 1s.
    - threshold: The threshold used for idealization.
    """
    
    # If no threshold is provided, estimate it
    if threshold is None:
        # Create a histogram to find the two main states (baseline noise vs. open channel current)
        hist, bin_edges = np.histogram(raw_data, bins=100)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Find peaks in the histogram (we expect two main peaks: closed and open)
        peaks, _ = find_peaks(hist, distance=10) # distance ensures peaks aren't too close
        
        if len(peaks) >= 2:
            # Sort peaks by amplitude to find the two highest (most populated states)
            sorted_peaks = sorted(peaks, key=lambda x: hist[x], reverse=True)[:2]
            # The threshold is roughly halfway between the two peaks
            val1 = bin_centers[sorted_peaks[0]]
            val2 = bin_centers[sorted_peaks[1]]
            threshold = (val1 + val2) / 2.0
            print(f"Auto-estimated threshold: {threshold:.3f}")
        else:
            # Fallback: simple mean if distinct peaks aren't found
            threshold = np.mean(raw_data)
            print(f"Warning: Could not distinct two states. Using mean as threshold: {threshold:.3f}")

    # Idealize the data based on the threshold
    # Assuming current is positive when open. If current is negative when open, invert the logic.
    # We will assume standard convention: open state has larger magnitude current.
    # Let's adjust logic based on the mean:
    
    mean_val = np.mean(raw_data)
    if mean_val > threshold:
         # Meaning the 'closed' state is higher (e.g. less negative current)
         idealized_states = (raw_data < threshold).astype(int)
    else:
         idealized_states = (raw_data > threshold).astype(int)


    if plot:
        # Plot a small segment (e.g., first 5000 points) to visually verify
        segment_len = min(5000, len(raw_data))
        t = np.arange(segment_len)
        
        fig, ax1 = plt.subplots(figsize=(10, 4))
        
        color = 'tab:blue'
        ax1.set_xlabel('Time (samples)')
        ax1.set_ylabel('Raw Current (pA)', color=color)
        ax1.plot(t, raw_data[:segment_len], color=color, alpha=0.6, label='Raw Data')
        ax1.axhline(threshold, color='red', linestyle='--', label=f'Threshold ({threshold:.2f})')
        ax1.tick_params(axis='y', labelcolor=color)
        
        ax2 = ax1.twinx()  
        color = 'tab:orange'
        ax2.set_ylabel('Idealized State (0/1)', color=color)  
        ax2.step(t, idealized_states[:segment_len], color=color, where='post', linewidth=2, label='Idealized')
        ax2.set_ylim(-0.1, 1.1)
        ax2.tick_params(axis='y', labelcolor=color)
        
        fig.tight_layout()  
        plt.title('Patch Clamp Idealization (First 5000 points)')
        plt.show()

    return idealized_states, threshold

# Example usage:
if __name__ == "__main__":
    # Simulate some raw patch clamp data (noisy baseline + noisy open states)
    np.random.seed(42)
    # Baseline noise (closed, mean=0 pA)
    baseline = np.random.normal(0, 0.5, 10000)
    # Open state events (mean=5 pA)
    events = np.zeros(10000)
    # Add random openings
    for _ in range(50):
        start = np.random.randint(0, 9500)
        length = np.random.randint(50, 300)
        events[start:start+length] = np.random.normal(5, 0.8, length)
        
    simulated_raw_data = baseline + events

    # 1. Idealize the data
    idealized_data, calc_threshold = idealize_patch_clamp(simulated_raw_data, plot=True)
    
    # 2. Save the idealized data for DFA analysis
    np.savetxt('idealized_sample.csv', idealized_data, fmt='%d')
    print("Idealized data saved to 'idealized_sample.csv'. You can now use this in dfa_hurst_analysis.py")
