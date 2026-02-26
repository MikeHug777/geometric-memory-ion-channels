# Geometric Memory in Ion Channels

**Pore Symmetry Determines Fractal Memory in Ion Channels**

Michael Hug (2026)

## The Formula

$$H = 1 - \frac{1}{B(n)}$$

where $B(n) = \frac{1}{n} \sum_{d \mid n} \varphi(n/d) \cdot 2^d$ is the Burnside orbit count for a binary ring of length $n$ under cyclic rotation $C_n$.

## Predictions

| Symmetry | B(n) | H predicted | Example channels |
|----------|------|-------------|------------------|
| C2 | 3 | 0.667 | TREK-2, TASK-3, Hv1 |
| C3 | 4 | 0.750 | ASIC1a, P2X, ENaC |
| C4 | 6 | 0.833 | BK, Kv, KcsA |
| C5 | 8 | 0.875 | nAChR, GABA_A, GlyR |
| C6 | 14 | 0.929 | Cx36, Orai/CRAC |

The formula retrodicts 30 years of published data from eight independent laboratories (1987-2024). TREK-2 (C2) matches to 1.0%, BK (C4) to 4.5%.

## Repository Structure

```
calculations/          41 Python scripts — all SI derivations and numerical verifications
data_analysis/         Patch-clamp analysis pipeline (DFA, idealization, batch processing)
figures/               Figure generation scripts for the paper
```

### Calculations

Each script in `calculations/` corresponds to one Supplementary Information section:

| Script | SI | Content |
|--------|----|---------|
| `burnside_orbit_topology.py` | SI-1 | Orbit enumeration, transition graphs, topology |
| `calculation_P_equipartition_derivation.py` | SI-2 | Core derivation: H = 1 - 1/B(n) |
| `calculation_AA_mckay_cartan_orbit.py` | SI-3 | McKay correspondence, Lie-algebraic foundation |
| `calculation_N_orbit_trap_model.py` | SI-4 | Bouchaud trap model on orbit graph |
| `calculation_S_necklace_return_time.py` | SI-5 | Necklace return-time distribution |
| `calculation_R_ergodic_orbit_equipartition.py` | SI-6 | Four ergodic justifications |
| `calculation_K2_DP_criticality_deep.py` | SI-7 | Directed percolation universality |
| `calculation_K_zeno_criticality.py` | SI-8 | Quantum Zeno and criticality |
| `calculation_J_XXZ_spin_ring.py` | SI-9 | XXZ spin chain, parameter-free H |
| `calculation_M_quantum_ring_spectrum.py` | SI-10 | Quantum ring spectrum |
| `calculation_O_quantum_temporal_correlator.py` | SI-11 | Quantum temporal correlator |
| `calculation_U_spectral_density_prediction.py` | SI-12 | Spectral density predictions |
| `calculation_Z_fano_factor_universal.py` | SI-13 | Universal Fano factor |
| `calculation_V_aging_ergodicity_breaking.py` | SI-14 | Aging and ergodicity breaking |
| `calculation_AC_coupled_burnside_traps.py` | SI-15 | Coupled network model |
| `calculation_W_multifractal_spectrum.py` | SI-16 | Multifractal spectrum |
| `calculation_AD_ising_ring_thermodynamics.py` | SI-17 | Ising ring thermodynamics |
| `calculation_T_temperature_phase_diagram.py` | SI-18 | Temperature phase diagram |
| `calculation_Q_information_thermodynamics.py` | SI-19 | Information thermodynamics |
| `calculation_AE_lowen_teich_retrodiction.py` | SI-20 | 30-year retrodiction |
| `calculation_E_thermodynamic_impossibility.py` | SI-21 | Classical impossibility |
| `calculation_G_H_I_breakthrough_attempts.py` | SI-22 | Spectral power excess |
| `calculation_L_ghost_of_symmetry.py` | SI-23 | Ghost of symmetry |
| `calculation_L2_simulate_H_from_kinetics.py` | SI-24 | H from kinetics |
| `calculation_X_anesthesia_prediction.py` | SI-25 | Anesthesia prediction |
| `calculation_Y_immune_system_prediction.py` | SI-26 | Immune system prediction |
| `calculation_F_autocorrelation_two_spin.py` | SI-27 | Two-spin autocorrelation |
| `calculation_AG_molien_partition.py` | SI-31 | Molien partition function |

### Data Analysis Pipeline

The `data_analysis/` directory contains a complete patch-clamp analysis pipeline:

- `idealize_patch_clamp.py` — Threshold-based idealization of raw recordings
- `dfa_hurst_analysis.py` — Detrended Fluctuation Analysis for Hurst exponent estimation
- `batch_patchclamp_analysis.py` — Batch processing across multiple recordings
- `run_pipeline_demo.py` — Demonstration of the full pipeline

## Requirements

- Python 3.8+
- NumPy
- SciPy
- Matplotlib (for figure scripts)

## Running

All calculation scripts are self-contained and produce terminal output:

```bash
python3 calculations/calculation_P_equipartition_derivation.py
```

## Supplementary Information

The full SI document (10 sections, derivations, numerical verifications, and data comparisons) is included in this repository:

- [`SI_geometric_memory.pdf`](SI_geometric_memory.pdf) — compiled PDF
- [`SI_geometric_memory.tex`](SI_geometric_memory.tex) — LaTeX source

## Paper

Hug, M. (2026). Pore Symmetry Determines Fractal Memory in Ion Channels. *bioRxiv* [preprint]. DOI: *forthcoming*

Pre-registration: [OSF](https://doi.org/10.17605/OSF.IO/BESWG)

## License

MIT License. See [LICENSE](LICENSE).
