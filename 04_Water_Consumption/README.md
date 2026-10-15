# Water Consumption — Fiber Photometry Analysis

This folder contains R scripts used to analyze mouth-aligned fiber photometry responses during water consumption. Recordings may also contain reaching events; the PETH extraction script selects mouth-related events for this analysis.

## Scripts

1. **`01_Water_Consumption_PETH.R`** — Reads recording metadata (`details.csv`), behavioral event files (`*_behavior.csv`), and photometry signals (`*_signal_df.csv`). Selects `body_part == "mouth"` in the metadata and behavioral events, applies a third-order Butterworth filter (`filtfilt`, normalized cutoff 0.1), and extracts ±5 s event-aligned traces. Each event is z-scored against its −4 to −1 s baseline. Events must be separated by at least 2 s from both neighboring events.
2. **`02_Water_Consumption_Graphs.R`** — Loads the combined PETH output, excludes `error` events, plots mouth-aligned responses, computes mean curves and trapezoidal AUCs per mouse and condition, and creates the figure combining PETHs and AUC comparisons. The supplied plotting code pools episodes for the PETH mean ± SEM, whereas AUCs are calculated from per-mouse mean curves. Paired t-tests compare time windows within each group and event type.

## Inputs

Set the project data path in each script before running. Expected inputs:

- `details.csv` with `video`, `mouse.ID`, `sex`, `group`, `phase`, and `body_part` columns.
- `*_behavior.csv` with `frame`, `event.type`, and an event label in `body_part` or `paw` (`mouth` for water consumption).
- `*_signal_df.csv` with time and photometry signal in the first two columns.

Both the metadata and the event files must identify mouth-related records.

## Run order

Run the scripts in numerical order:

```r
source("01_Water_Consumption_PETH.R")
source("02_Water_Consumption_Graphs.R")
```

Each script clears the R workspace; run them individually after configuring their data paths.

## Reproducibility notes

- Photometry and behavioral frame periods are both configured as `1/60` s.
- PETH windows extend from −5 to +5 s around each event.
- The graph script reports paired within-group comparisons between AUC windows.
