# ============================================================
# Reaching Test — Behavioral Performance Analysis
# ============================================================
#
# Description:
#   Quantify reaching performance across experimental groups and test days using manually annotated behavioral events.
#
#   - Select Day 1 and Day 4 recordings
#   - Determine preferred paw and reaching outcomes
#   - Quantify event counts, outcome ratios, and target distances
#   - Export per-recording summaries and individual events
#   - Reuse for other experiments by changing EXPERIMENT
#
# Input:
#   <DATA_DIR>/<EXPERIMENT>/reaching_test_details_<EXPERIMENT>.csv
#   <DATA_DIR>/<EXPERIMENT>/**/*_behavior.csv
#
# Output:
#   <DATA_DIR>/<EXPERIMENT>/results_mean_<EXPERIMENT>.csv
#   <DATA_DIR>/<EXPERIMENT>/results_all_<EXPERIMENT>.csv
#
# Author: Alessandro Contestabile
# ============================================================

# ============================================================
# 1. SETUP
# ============================================================
rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

# ============================================================
# 2. PATHS AND PARAMETERS
# ============================================================
DATA_DIR <- "path/to/reaching_test"
EXPERIMENT <- "anxa1_tet"  # Change this to analyze another experiment
EXP_DIR <- file.path(DATA_DIR, EXPERIMENT)
DETAILS_FILE <- file.path(DATA_DIR, "details.csv")

FRAME_DURATION <- 1 / 60
PHASES <- c("test1", "test4")
GROUPS <- c("ctrl", "tet")  # Adapt to the groups in the selected experiment
GROUP_COLORS <- c("ctrl" = "grey", "tet" = "#ffb437")  # For downstream plotting
TARGET_REFERENCE_CM <- 2
DISTANCE_THRESHOLD_CM <- 0.25

stopifnot(dir.exists(EXP_DIR), file.exists(DETAILS_FILE))

# ============================================================
# 3. LOAD METADATA AND DISCOVER FILES
# ============================================================
details <- readr::read_csv(DETAILS_FILE, show_col_types = FALSE) %>%
  dplyr::filter(
    experiment == EXPERIMENT,
    phase %in% PHASES
  )

required_cols <- c("video", "mouse.ID", "day", "group", "phase")
stopifnot(all(required_cols %in% names(details)))

behavior_files <- list.files(EXP_DIR, pattern = "_behavior\\.csv$", recursive = TRUE, full.names = TRUE)
if (!length(behavior_files)) stop("No *_behavior.csv files found.")

# ============================================================
# 4. PROCESS INDIVIDUAL RECORDINGS
# ============================================================
xy_dataset <- list()

summary_cols <- c(
  "prefered.paw", "prefered_tot", "prefered_perc",
  "rewarded_tot", "rewarded_perc", "rewarded_dist",
  "errors_tot", "errors_perc", "errors_dist",
  "empty_tot", "empty_perc", "empty_dist",
  "mean_distance", "perc_thr"
)
for (nm in summary_cols) details[[nm]] <- if (nm == "prefered.paw") NA_character_ else NA_real_

for (behavior_file in behavior_files) {
  # The recording identifier is the behavior filename without its suffix.
  video <- sub("_behavior\\.csv$", "", basename(behavior_file))
  info_idx <- which(as.character(details$video) == video)
  if (!length(info_idx)) next
  if (length(info_idx) > 1L) {
    warning("Multiple metadata rows for video ", video, "; using the first.")
    info_idx <- info_idx[1]
  }
  info <- details[info_idx, , drop = FALSE]

  df <- read.csv(behavior_file, stringsAsFactors = FALSE)
  if (ncol(df) < 10) {
    warning("Skipping ", video, ": fewer than 10 behavioral columns.")
    next
  }
  df <- df[, 1:10]
  names(df) <- c("frame", "paw", "event.type", "target.x", "target.y",
                 "paw.x", "paw.y", "reference.px.", "distance.px.", "distance.cm.")

  # Preserve original complete-case filtering and paw selection.
  df <- na.omit(df)
  df <- df[df$paw %in% c("left", "right", "L", "R"), , drop = FALSE]
  events_tot <- nrow(df)
  if (!events_tot) {
    warning("Skipping ", video, ": no valid paw events.")
    next
  }

  left_paw <- sum(df$paw %in% c("left", "L")) / events_tot
  right_paw <- sum(df$paw %in% c("right", "R")) / events_tot
  prefered_paw <- ifelse(left_paw > right_paw, "left", "right")

  df <- df %>%
    mutate(
      event.type = case_when(
        event.type %in% c("correct", "semi-drop", "C") ~ "Rewarded",
        event.type %in% c("drop", "missed", "miss", "M") ~ "Errors",
        event.type %in% c("empty", "E") ~ "Empty",
        TRUE ~ as.character(event.type)
      ),
      paw = case_when(
        paw %in% c("L", "Left", "left") ~ "left",
        paw %in% c("R", "Right", "right") ~ "right",
        TRUE ~ as.character(paw)
      ),
      nor.paw.x = ((paw.x - target.x) * TARGET_REFERENCE_CM) / reference.px.,
      nor.paw.y = ((paw.y - target.y) * TARGET_REFERENCE_CM) / reference.px.,
      mouse.ID = as.character(info$mouse.ID[1]),
      video = video,
      day = info$day[1],
      group = info$group[1],
      phase = info$phase[1]
    )

  prefered_tot <- sum(df$paw == prefered_paw)
  prefered_perc <- prefered_tot / events_tot

  # Original analysis: outcome counts include all paw events, while
  # outcome percentages use preferred-paw event count as denominator.
  # This is retained for reproduction; ratios may exceed 1.
  outcome_summary <- function(outcome) {
    selected <- df$event.type == outcome
    c(total = sum(selected),
      percentage = if (prefered_tot > 0) sum(selected) / prefered_tot else NA_real_,
      distance = mean(df$distance.cm.[selected], na.rm = TRUE))
  }
  rewarded <- outcome_summary("Rewarded")
  errors <- outcome_summary("Errors")
  empty <- outcome_summary("Empty")

  values <- list(
    prefered.paw = prefered_paw,
    prefered_tot = prefered_tot,
    prefered_perc = prefered_perc,
    rewarded_tot = rewarded["total"],
    rewarded_perc = rewarded["percentage"],
    rewarded_dist = rewarded["distance"],
    errors_tot = errors["total"],
    errors_perc = errors["percentage"],
    errors_dist = errors["distance"],
    empty_tot = empty["total"],
    empty_perc = empty["percentage"],
    empty_dist = empty["distance"],
    mean_distance = mean(df$distance.cm., na.rm = TRUE),
    perc_thr = mean(df$distance.cm. < DISTANCE_THRESHOLD_CM, na.rm = TRUE) * 100
  )
  for (nm in names(values)) details[info_idx, nm] <- unname(values[[nm]])

  xy_dataset[[length(xy_dataset) + 1L]] <- df
  message("Processed: ", video)
}

# ============================================================
# 5. PREPARE SUMMARY DATA
# ============================================================
details <- details %>%
  mutate(
    mouse.ID = as.character(mouse.ID),
    phase = paste("Day", sub("^test", "", phase)),
    group = factor(group, levels = GROUPS),
    phase = factor(phase, levels = c("Day 1", "Day 4"))
  ) %>%
  filter(!is.na(phase)) %>%
  mutate(across(ends_with("_tot") | ends_with("_perc"), ~ replace_na(.x, 0)))

# ============================================================
# 6. EXPORT RESULTS
# ============================================================
write.csv(details, file.path(EXP_DIR, paste0("results_mean_", EXPERIMENT, ".csv")), row.names = FALSE)

if (!length(xy_dataset)) stop("No behavioral events processed; summary saved but event output is empty.")
all_events <- bind_rows(xy_dataset)
write.csv(all_events, file.path(EXP_DIR, paste0("results_all_", EXPERIMENT, ".csv")), row.names = FALSE)

message("Reaching test behavioral analysis completed.")
