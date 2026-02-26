"""
Data Conversion — Wawrzkiewicz-Jalowiecka BK Channel Dataset (Zenodo 7817601)
================================================================================
Author:  Michael Hug
Date:    2026-02-22
Project: Interface First — Geometric Memory
Paper:   Paper 2, Analysis Pipeline (Data)

Original contribution:
  Converts the Wawrzkiewicz-Jalowiecka et al. (2024) BK single-channel
  dataset (Zenodo DOI: 10.5281/zenodo.7817601) from Excel format to
  CSV time series for DFA analysis. This dataset provides the C4 (BK)
  and C2 (TREK-2) recordings used for retrodiction validation.

Dependencies: numpy, pandas
"""
import csv
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    base = Path("/Users/mikehug/Documents/personality-overview-core/pcm_model/data/Forschung/Interface_First/calculations")
    xlsx_dir = base / "external_data" / "extracted" / "Data"
    real_data_dir = base / "real_data"
    real_data_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = base / "manifest_real_run_zenodo.csv"

    xlsx_files = sorted(xlsx_dir.glob("*_raw.xlsx"))
    if not xlsx_files:
        raise FileNotFoundError(f"No *_raw.xlsx files found in {xlsx_dir}")

    rows = []
    for xlsx_path in xlsx_files:
        # Read without header because these files store numeric data directly.
        df = pd.read_excel(xlsx_path, header=None)

        if df.shape[1] < 2:
            # fallback: only one column available
            current = pd.to_numeric(df.iloc[:, 0], errors="coerce").dropna().to_numpy(dtype=float)
            chosen_col = 0
        else:
            # Assume col 0 is time and cols 1..N are current traces.
            current_block = df.iloc[:, 1:]
            current_block = current_block.apply(pd.to_numeric, errors="coerce")

            # Choose the trace with highest variance (typically richest gating dynamics).
            stds = current_block.std(axis=0, skipna=True).to_numpy(dtype=float)
            best_rel_idx = int(np.nanargmax(stds))
            chosen_col = best_rel_idx + 1

            current = current_block.iloc[:, best_rel_idx].dropna().to_numpy(dtype=float)

        out_name = f"{xlsx_path.stem}_col{chosen_col:02d}.csv"
        out_path = real_data_dir / out_name
        np.savetxt(out_path, current, fmt="%.8f")

        rows.append(
            {
                "filename": out_name,
                "channel_name": xlsx_path.stem,
                "symmetry": "C4",
                "expected_h": "0.833",
                "threshold": "",
            }
        )

    with manifest_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["filename", "channel_name", "symmetry", "expected_h", "threshold"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Converted {len(rows)} files into {real_data_dir}")
    print(f"Manifest written to: {manifest_path}")


if __name__ == "__main__":
    main()
