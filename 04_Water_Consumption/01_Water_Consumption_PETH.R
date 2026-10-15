# ============================================================
# Water Consumption — Fiber Photometry PETH
# ============================================================
#
# Description:
#   Align fiber photometry signals to mouth-related water consumption events in recordings that may also contain reaching behavior.
#
#   - Select mouth recordings and events from shared datasets
#   - Exclude events with neighboring events within 2 seconds
#   - Apply third-order, zero-phase Butterworth filtering
#   - Extract +/-5-second peri-event photometry signals
#   - Compute event-wise z-scores using a -4 to -1 s baseline
#   - Export per-video PETHs and combined analysis tables
#
# Input:
#   details.csv (video, mouse.ID, sex, group, phase, body_part)
#   <video>_signal_df.csv
#   <video>_behavior.csv (frame, event.type, body_part or paw)
#
# Output:
#   results/<mouse.ID>/<video>_PETH.csv
#   results/PETH_all.csv
#   results/PETH_all.rds
#   results/run_summary.csv
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
  library(stringr)
  library(tibble)
  library(signal)
})

# ============================================================
# 2. PATHS AND PARAMETERS
# ============================================================

DATA_DIR <- "path/to/data"
OUT_DIR <- file.path(DATA_DIR, "results")
DETAILS_FILE <- file.path(DATA_DIR, "details.csv")

FRAME_DF <- 1 / 60
FRAME_BEHAVIOR <- 1 / 60
PETH_WINDOW <- 5
BASELINE_WINDOW <- c(-4, -1)
MIN_EVENT_GAP <- 2
FILTER_ORDER <- 3
FILTER_CUTOFF <- 0.1

stopifnot(dir.exists(DATA_DIR), file.exists(DETAILS_FILE))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 3. LOAD AND FILTER METADATA
# ============================================================

normalize_video <- function(x) {
  x <- tools::file_path_sans_ext(basename(as.character(x)))
  x <- sub("_signal_df$", "", x)
  gsub(":", "_", x)
}

details <- read_csv(DETAILS_FILE, show_col_types = FALSE)
required_cols <- c("video", "mouse.ID", "sex", "group", "phase", "body_part")
stopifnot(all(required_cols %in% names(details)))

details <- details %>%
  mutate(video = normalize_video(video), body_part = str_to_lower(str_trim(as.character(body_part)))) %>%
  filter(body_part == "mouth")

if (!nrow(details)) stop("No mouth entries in details.csv (body_part == 'mouth').")


# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================

# Keep only events isolated by at least MIN_EVENT_GAP on both sides.
filter_isolated_events <- function(df, min_gap = 2) {
  if (!nrow(df)) return(df)
  df <- df[order(df$time), , drop = FALSE]
  gaps <- diff(df$time)
  df[c(Inf, gaps) >= min_gap & c(gaps, Inf) >= min_gap, , drop = FALSE]
}

# Event-wise z-score using the specified pre-event baseline.
zscore_window <- function(signal_vec, time_vec, event_time, baseline_window) {
  idx <- which(time_vec >= event_time + baseline_window[1] &
               time_vec <= event_time + baseline_window[2])
  if (length(idx) < 10) {
    mu <- mean(signal_vec, na.rm = TRUE)
    sdv <- sd(signal_vec, na.rm = TRUE)
  } else {
    mu <- mean(signal_vec[idx], na.rm = TRUE)
    sdv <- sd(signal_vec[idx], na.rm = TRUE)
  }
  if (!is.finite(sdv) || sdv == 0) sdv <- 1e-9
  (signal_vec - mu) / sdv
}

# ============================================================
# 5. PROCESS RECORDINGS
# ============================================================

signal_files <- list.files(DATA_DIR, pattern = "_signal_df\\.csv$", recursive = TRUE, full.names = TRUE)
signal_files <- signal_files[!startsWith(normalizePath(signal_files, mustWork = FALSE),
                                        normalizePath(OUT_DIR, mustWork = FALSE))]
if (!length(signal_files)) stop("No *_signal_df.csv files found.")

all_peth <- list()
run_summary <- list()

for (signal_file in signal_files) {
  video <- normalize_video(signal_file)
  meta <- details %>% filter(.data$video == video)
  if (!nrow(meta)) next
  if (nrow(meta) > 1) {
    warning("Multiple mouth metadata rows for ", video, "; using the first.")
    meta <- slice_head(meta, n = 1)
  }

  # ------------------------------------------------------------
  # 5.1 LOAD AND FILTER BEHAVIORAL EVENTS
  # ------------------------------------------------------------

  behavior_file <- sub("_signal_df\\.csv$", "_behavior.csv", signal_file)
  if (!file.exists(behavior_file)) {
    warning("Missing behavior file for ", video)
    next
  }

  behavior <- read_csv(behavior_file, show_col_types = FALSE)
  if (!all(c("frame", "event.type") %in% names(behavior))) {
    warning("Missing frame/event.type in ", video)
    next
  }

  # Mixed recordings may include reaching and mouth events.
  if ("body_part" %in% names(behavior)) {
    behavior <- behavior %>% filter(str_to_lower(str_trim(as.character(body_part))) == "mouth")
  } else if ("paw" %in% names(behavior)) {
    behavior <- behavior %>% filter(str_to_lower(str_trim(as.character(paw))) == "mouth")
  } else {
    warning("Cannot identify mouth events in ", video, "; skipping mixed recording.")
    next
  }

  behavior <- behavior %>%
    mutate(frame = suppressWarnings(as.numeric(frame)), time = frame * FRAME_BEHAVIOR) %>%
    filter(is.finite(time)) %>% arrange(time)

  n_before <- nrow(behavior)
  behavior <- filter_isolated_events(behavior, MIN_EVENT_GAP)

  # ------------------------------------------------------------
  # 5.2 LOAD AND FILTER PHOTOMETRY
  # ------------------------------------------------------------

  signal_df <- read_csv(signal_file, show_col_types = FALSE)
  if (ncol(signal_df) < 2) next
  names(signal_df)[1:2] <- c("time", "dF")
  signal_df <- signal_df %>%
    transmute(time = suppressWarnings(as.numeric(time)), dF = suppressWarnings(as.numeric(dF))) %>%
    filter(is.finite(time), is.finite(dF)) %>%
    arrange(time) %>% distinct(time, .keep_all = TRUE)
  if (nrow(signal_df) < 5) next

  bf <- signal::butter(FILTER_ORDER, FILTER_CUTOFF)
  signal_df$dF_filt <- tryCatch(as.numeric(signal::filtfilt(bf, signal_df$dF)),
    error = function(e) { warning("Filter failed for ", video, ": ", e$message); signal_df$dF })

  # ------------------------------------------------------------
  # 5.3 EXTRACT AND Z-SCORE PERI-EVENT SIGNALS
  # ------------------------------------------------------------

  halfwin <- round(PETH_WINDOW / FRAME_DF)
  nor_time <- seq(-halfwin, halfwin) * FRAME_DF
  peth_list <- vector("list", nrow(behavior))

  for (i in seq_len(nrow(behavior))) {
    onset <- behavior$time[i]
    if (onset <= min(signal_df$time) + PETH_WINDOW ||
        onset >= max(signal_df$time) - PETH_WINDOW) next

    center <- which.min(abs(signal_df$time - onset))
    idx <- (center - halfwin):(center + halfwin)
    if (min(idx) < 1 || max(idx) > nrow(signal_df)) next
    segment <- signal_df[idx, ]
    if (any(!is.finite(segment$dF_filt))) next

    peth_list[[i]] <- tibble(
      nor.time = nor_time,
      zscore = zscore_window(segment$dF_filt, segment$time, onset, BASELINE_WINDOW),
      event.type = as.character(behavior$event.type[i]),
      body_part = "mouth",
      episode = i
    )
  }

  # ------------------------------------------------------------
  # 5.4 ADD METADATA AND SAVE VIDEO-LEVEL RESULTS
  # ------------------------------------------------------------

  peth <- bind_rows(peth_list)
  run_summary[[length(run_summary) + 1]] <- tibble(
    video = video, mouse.ID = as.character(meta$mouse.ID[1]), group = as.character(meta$group[1]),
    n_mouth_events = n_before, n_isolated = nrow(behavior),
    n_kept = length(Filter(Negate(is.null), peth_list))
  )
  if (!nrow(peth)) next

  peth <- peth %>% mutate(
    mouse.ID = as.character(meta$mouse.ID[1]), video = video,
    phase = as.character(meta$phase[1]), group = as.character(meta$group[1]),
    sex = as.character(meta$sex[1])
  )

  mouse_dir <- file.path(OUT_DIR, as.character(meta$mouse.ID[1]))
  dir.create(mouse_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(peth, file.path(mouse_dir, paste0(video, "_PETH.csv")))
  all_peth[[length(all_peth) + 1]] <- peth
  message("Processed ", video, ": ", run_summary[[length(run_summary)]]$n_kept, " events")
}


# ============================================================
# 6. COMBINE AND EXPORT
# ============================================================

if (!length(all_peth)) stop("No mouth PETH events were generated.")

PETH_all <- bind_rows(all_peth) %>% mutate(
  event.type = case_when(
    str_to_lower(event.type) %in% c("c", "correct", "semi-drop", "semidrop") ~ "rewarded",
    str_to_lower(event.type) %in% c("drop", "missed", "miss", "e", "m") ~ "error",
    str_to_lower(event.type) %in% c("empty", "emopty") ~ "empty",
    TRUE ~ event.type
  )
)

write_csv(PETH_all, file.path(OUT_DIR, "PETH_all.csv"))
saveRDS(PETH_all, file.path(OUT_DIR, "PETH_all.rds"))
write_csv(bind_rows(run_summary), file.path(OUT_DIR, "run_summary.csv"))
message("Water consumption PETH analysis completed.")