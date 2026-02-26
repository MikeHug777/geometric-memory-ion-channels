# Batch Patch-Clamp Workflow

This workflow runs **idealization + DFA** for multiple raw patch-clamp CSV files.

## 1) Input format
- Put raw current files (`*.csv`) in one folder.
- Each file should contain one numeric column (current trace).

Optional: use `manifest_template.csv` to set channel metadata.

## 2) Run
```bash
/opt/homebrew/bin/python3 batch_patchclamp_analysis.py \
  --input-dir /ABS/PATH/to/raw_data \
  --output-dir /ABS/PATH/to/output \
  --manifest /ABS/PATH/to/manifest.csv
```

Without manifest:
```bash
/opt/homebrew/bin/python3 batch_patchclamp_analysis.py \
  --input-dir /ABS/PATH/to/raw_data \
  --output-dir /ABS/PATH/to/output
```

## 3) Outputs
For each input file:
- `*_idealized.csv` (0/1 states)
- `*_dfa.png` (log-log DFA fit)

Global:
- `summary.csv` with:
  - filename
  - channel_name
  - symmetry
  - threshold
  - expected_h
  - calculated_h
  - deviation_percent

## 4) Notes
- Expected H defaults from symmetry if not manually provided:
  - C2: 0.667
  - C3: 0.750
  - C4: 0.833
  - C5: 0.875
  - C6: 0.929
- If threshold is omitted, it is auto-estimated from histogram peaks.
