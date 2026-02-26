"""
Batch Patch-Clamp Analysis with Cn Symmetry Inference
=======================================================
Author:  Michael Hug
Date:    2026-02-22
Project: Interface First — Geometric Memory
Paper:   Paper 2, Analysis Pipeline

Original contribution:
  Batch processing pipeline that automatically infers Cn symmetry class
  from channel name, computes DFA Hurst exponents, and compares measured
  H values against Burnside predictions H = 1 - 1/B(n). Enables
  systematic validation across multiple recordings and channel types.

Dependencies: numpy, matplotlib, argparse
"""
import argparse
import csv
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from idealize_patch_clamp import idealize_patch_clamp
from dfa_hurst_analysis import dfa


def expected_h_from_symmetry(symmetry: str):
    """
    Burnside-based expected H values used in your framework.
    """
    lookup = {
        "C2": 2.0 / 3.0,
        "C3": 0.750,
        "C4": 5.0 / 6.0,
        "C5": 0.875,
        "C6": 13.0 / 14.0,
    }
    return lookup.get(symmetry)


def infer_symmetry_from_name(name: str):
    upper = name.upper()
    for tag in ["C2", "C3", "C4", "C5", "C6"]:
        if tag in upper:
            return tag
    return ""


def load_manifest(manifest_path: Path):
    """
    Optional manifest CSV columns:
    filename,channel_name,symmetry,expected_h,threshold
    """
    if not manifest_path.exists():
        return {}

    result = {}
    with manifest_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            filename = row.get("filename", "").strip()
            if not filename:
                continue
            result[filename] = {
                "channel_name": row.get("channel_name", "").strip(),
                "symmetry": row.get("symmetry", "").strip().upper(),
                "expected_h": row.get("expected_h", "").strip(),
                "threshold": row.get("threshold", "").strip(),
            }
    return result


def save_dfa_plot(scales, fluctuations, h_value, expected_h, title, out_path: Path):
    plt.figure(figsize=(8, 6))
    plt.loglog(scales, fluctuations, "o", label="DFA fluctuation")
    fit = np.polyfit(np.log10(scales), np.log10(fluctuations), 1)
    fit_line = 10 ** np.polyval(fit, np.log10(scales))
    plt.loglog(scales, fit_line, "-", label=f"Fit H={h_value:.3f}")

    if expected_h is not None:
        plt.title(f"{title} | expected H={expected_h:.3f}")
    else:
        plt.title(title)

    plt.xlabel("Scale")
    plt.ylabel("Fluctuation")
    plt.grid(True, which="both", ls="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def parse_float_or_none(value):
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    value = str(value).strip()
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def process_file(csv_path: Path, metadata: dict, output_dir: Path, already_idealized: bool = False):
    raw = np.loadtxt(csv_path)
    if raw.ndim > 1:
        raw = raw[:, 0]

    channel_name = metadata.get("channel_name") or csv_path.stem
    symmetry = (metadata.get("symmetry") or infer_symmetry_from_name(csv_path.stem)).upper()

    expected_h = parse_float_or_none(metadata.get("expected_h"))
    if expected_h is None and symmetry:
        expected_h = expected_h_from_symmetry(symmetry)

    threshold = parse_float_or_none(metadata.get("threshold"))

    if already_idealized:
        # Treat input as state trace; binarize multi-level idealization (e.g. 0/1/2) to open-vs-closed.
        idealized = (raw > 0).astype(int)
        used_threshold = None
    else:
        idealized, used_threshold = idealize_patch_clamp(raw, threshold=threshold, plot=False)
    scales, fluctuations, h_value = dfa(idealized)

    error_pct = None
    if expected_h is not None and expected_h != 0:
        error_pct = abs(h_value - expected_h) / expected_h * 100.0

    idealized_out = output_dir / f"{csv_path.stem}_idealized.csv"
    np.savetxt(idealized_out, idealized.astype(int), fmt="%d")

    plot_out = output_dir / f"{csv_path.stem}_dfa.png"
    save_dfa_plot(
        scales=scales,
        fluctuations=fluctuations,
        h_value=h_value,
        expected_h=expected_h,
        title=channel_name,
        out_path=plot_out,
    )

    return {
        "filename": csv_path.name,
        "channel_name": channel_name,
        "symmetry": symmetry,
        "threshold": used_threshold,
        "expected_h": expected_h,
        "calculated_h": h_value,
        "deviation_percent": error_pct,
        "idealized_file": idealized_out.name,
        "dfa_plot": plot_out.name,
    }


def main():
    parser = argparse.ArgumentParser(description="Batch idealization + DFA for patch-clamp CSV files")
    parser.add_argument("--input-dir", required=True, help="Directory with raw CSV files")
    parser.add_argument("--output-dir", required=True, help="Directory for outputs")
    parser.add_argument(
        "--manifest",
        default="",
        help="Optional manifest CSV (filename,channel_name,symmetry,expected_h,threshold)",
    )
    parser.add_argument(
        "--glob",
        default="*.csv",
        help="Input glob pattern (default: *.csv)",
    )
    parser.add_argument(
        "--already-idealized",
        action="store_true",
        help="Skip threshold idealization and treat input traces as idealized states.",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest_map = {}
    if args.manifest:
        manifest_map = load_manifest(Path(args.manifest).expanduser().resolve())

    files = sorted(input_dir.glob(args.glob))
    if not files:
        raise FileNotFoundError(f"No files found in {input_dir} with glob '{args.glob}'")

    rows = []
    for file_path in files:
        # Skip already idealized files by default
        if file_path.stem.endswith("_idealized"):
            continue
        metadata = manifest_map.get(file_path.name, {})
        try:
            row = process_file(file_path, metadata, output_dir, already_idealized=args.already_idealized)
            rows.append(row)
            print(f"[OK] {file_path.name}: H={row['calculated_h']:.3f}")
        except Exception as exc:
            print(f"[ERROR] {file_path.name}: {exc}")

    summary_path = output_dir / "summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "filename",
                "channel_name",
                "symmetry",
                "threshold",
                "expected_h",
                "calculated_h",
                "deviation_percent",
                "idealized_file",
                "dfa_plot",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nDone. Summary written to: {summary_path}")


if __name__ == "__main__":
    main()
