# Reaching Test — Behavioral Analysis

This folder contains two complementary scripts for quantifying mouse reaching behavior from manually annotated events and DeepLabCut tracking data.

## Scripts

- **`01_Reaching_Test_Behavior.R`** — Reads `details.csv`, selects the experiment specified by `EXPERIMENT` and the phases in `PHASES` (default: `test1` and `test4`), and processes annotated `*_behavior.csv` files. Computes preferred-paw usage, rewarded/error/empty event counts and ratios, and distances from the target. Exports `results_mean_<EXPERIMENT>.csv` and `results_all_<EXPERIMENT>.csv`.
- **`02_Reaching_Test_DLC.py`** — Processes DeepLabCut CSV files containing snout, palm, and spout coordinates (mouth coordinates are optional). Extracts movement bouts, palm advance/grasp/retract metrics, trajectory and grasp measures, and generates per-recording PDF plots. Exports `trace_distance_summary.csv` and per-bout phase tables.

## Requirements

- **R:** `dplyr`, `readr`, `tidyr`, `stringr`
- **Python:** `pandas`, `numpy`, `matplotlib`

## Usage

1. In the R script, set `DATA_DIR`, `EXPERIMENT`, `PHASES`, and `GROUPS`. The metadata file is `<DATA_DIR>/details.csv`; annotated behavior files are searched under `<DATA_DIR>/<EXPERIMENT>/`.
2. In the Python script, set `DATA_DIR`, `CSV_PATTERN`, and `OUT_DIR` to the relevant DeepLabCut CSV directory and output location.
3. Run each script independently. To reuse the R analysis for another experimental cohort, change `EXPERIMENT` and the group settings as needed.

**Analysis notes:** The R script preserves the original outcome-ratio definition (outcomes counted across paw events, divided by preferred-paw event count). The Python script uses pixel-based thresholds, including `ACCURACY_RADIUS_PX`; review these settings against the imaging scale before interpreting spatial accuracy. Both scripts expect the input files to follow the formats described in their headers and code.
