# Open Field Test

R scripts used for open-field behavioral analyses, VAME motif characterization, fiber-photometry alignment, and figure generation.

## Setup

Install the R packages listed at the beginning of each script. Some preprocessing scripts additionally require a working Python environment (accessed through `reticulate`) for VAME/DLC files.

Before running, update the input/output paths in the **Settings** or **Paths** section of each script (typically `path/to/VAME_project` or `path/to/VAME_project/results`). Input recordings, annotations, and intermediate data files are not included here.

## Recommended execution order

1. **`00_VAME_Motif_Selection.R`** — Inspect motif frequencies and select motifs for analysis.
2. **`01_VAME_Dendrogram.R`** — Examine motif similarity and the motif grouping used downstream.
3. **`02_Data_Alignment.R`** — Align VAME motifs, behavioral tracking, and photometry; export `*_aligned.csv` and recording summaries.
4. **`03_PETH_Preprocessing.R`** — Generate motif/event-aligned photometry summaries.
5. **`Figure_5D-H.R`** — Generate Figure 5 analyses and the motif-level intermediate tables used by subsequent scripts.
6. **`Figure_6C-I.R`** — Generate Figure 6 behavioral summaries using the preprocessing outputs and motif-speed table.
7. **`00_VAME_Feature_Extraction.R`** — Build `vame_sequence_features.csv` from aligned recordings and the motif-speed table.
8. **`01_XGBoost_Classification_Figure_6J-K.R`** — Train the TeTxLC classifier and evaluate the other experimental cohorts (Figure 6J–K).
9. **`Figure_6M.R`** — Plot behavioral transition networks using the extracted sequence features and `summary_annotations.csv`.

## Extended Data figures

- **`Extended_Data_Figure_6A.R`** — Photometry signal versus speed, with circular-shuffle controls; requires aligned recordings.
- **`Extended_Data_Figure_6B.R`** — Example aligned photometry, kinematics, and motif ethogram; requires aligned recordings.
- **`Extended_Data_Figure_6CDE.R`** — Motif episode durations and short-episode fractions; requires `motif_episodes.csv` and `mean_speed_per_motif_overall.csv` from `Figure_5D-H.R`.
- **`Extended_Data_Figure_8I.R`** — Speed distributions across experimental cohorts; requires aligned recordings.

Generated figures and summary tables are written to the configured `results/analysis` directory unless otherwise specified in a script.
