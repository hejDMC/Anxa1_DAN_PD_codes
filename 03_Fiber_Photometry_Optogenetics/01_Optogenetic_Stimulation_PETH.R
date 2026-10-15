# ============================================================
# Optogenetic Stimulation — Fiber Photometry PETH
# ============================================================
#
# Description:
#   Align fiber photometry signals to optogenetic stimulation onset during open-field recordings.
#
#   - Identify stimulation onsets from Bonsai TTL recordings
#   - Detect stimulation bursts and pulse positions
#   - Extract peri-stimulation photometry signals
#   - Normalize fluorescence to recording maximum
#   - Plot mean +/- SEM for each pulse position
#   - Export PETH data and figures per recording
#
# Input:
#   oft_fp_opto_details.csv
#   <experiment>/videos/<video>.csv
#   <experiment>/signal/<video>_signal_df.csv
#
# Output:
#   results/<video>/<video>_PETH.csv
#   results/<video>/<video>_PETH.pdf
#
# Author: Alessandro Contestabile
# ============================================================

# ============================================================
# 1. SETUP
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(signal)
})


# ============================================================
# 2. PATHS AND PARAMETERS
# ============================================================

DATA_DIR <- "path/to/oft/opto+dLight"
OUT_DIR <- file.path(DATA_DIR, "results")

DETAILS_FILE <- file.path(DATA_DIR, "oft_fp_opto_details.csv")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FRAME_BEHAVIOR <- 1 / 30
FRAME_DF <- 1 / 60

PETH_WINDOW <- 2
PROTOCOL_ANALYSIS <- "1 sec, 3 sec"          # Change here if you want another protocol in details.csv

MIN_RECORDING_TIME <- 60
BURST_INTERVAL <- 4

FILTER_ORDER <- 3
FILTER_CUTOFF <- 0.1

BIN_WIDTH <- 0.2
SKIP_EXISTING <- TRUE

STIM_COLOR <- "steelblue"

# ============================================================
# 3. LOAD EXPERIMENTAL METADATA
# ============================================================

details <- read.csv(DETAILS_FILE, stringsAsFactors = FALSE)

required_cols <- c("mouse.ID", "group", "video", "experiment", "start.fp", "protocol", "intenisty", "frequency")
stopifnot(all(required_cols %in% colnames(details)))

details <- details %>%
  filter(protocol == PROTOCOL_ANALYSIS, !is.na(start.fp), !is.na(video))

cat("Selected recordings:", nrow(details), "\n")


# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================

# Calculate mean and SEM for each time bin and pulse position
summarize_peth <- function(df) {
  df %>%
    mutate(time_bin = round(nor.time / BIN_WIDTH) * BIN_WIDTH) %>%
    group_by(episode.number, time_bin) %>%
    summarise(
      mean_signal = mean(normalized_dF, na.rm = TRUE),
      sd_signal = sd(normalized_dF, na.rm = TRUE),
      n = sum(!is.na(normalized_dF)),
      sem = ifelse(n > 1, sd_signal / sqrt(n), NA_real_),
      .groups = "drop"
    )
}


# ============================================================
# 5. PROCESS EACH RECORDING
# ============================================================

for (f in seq_len(nrow(details))) {
  
  info <- details[f, ]
  
  video <- info$video
  experiment <- info$experiment
  
  video_file <- file.path(DATA_DIR, experiment, "videos", paste0(video, ".csv"))
  signal_file <- file.path(DATA_DIR, experiment, "signal", paste0(video, "_signal_df.csv"))
  
  video_out <- file.path(OUT_DIR, video)
  dir.create(video_out, recursive = TRUE, showWarnings = FALSE)
  
  csv_out <- file.path(video_out, paste0(video, "_PETH.csv"))
  pdf_out <- file.path(video_out, paste0(video, "_PETH.pdf"))
  
  if (SKIP_EXISTING && file.exists(csv_out) && file.exists(pdf_out)) {
    message("Skipping existing: ", video)
    next
  }
  
  if (!file.exists(video_file) || !file.exists(signal_file)) {
    warning("Missing input files for: ", video)
    next
  }
  
  message("Processing: ", video)
  
  # ============================================================
  # 5.1 LOAD BONSAI AND PHOTOMETRY DATA
  # ============================================================
  
  bonsai <- read.csv(video_file)
  dF <- read.csv(signal_file)
  
  if (ncol(bonsai) != 3 || ncol(dF) != 2) {
    warning("Unexpected input format: ", video)
    next
  }
  
  colnames(bonsai) <- c("timestamp", "time", "in")
  colnames(dF) <- c("time", "dF")
  
  bonsai$time <- as.numeric(bonsai$time) * FRAME_BEHAVIOR
  bonsai$in <- tolower(as.character(bonsai$in)) == "true"
  
  dF$time <- as.numeric(dF$time)
  dF$dF <- as.numeric(dF$dF)
  
  bonsai <- bonsai %>%
    filter(time > as.numeric(info$start.fp) * FRAME_BEHAVIOR + MIN_RECORDING_TIME)
  
  if (nrow(bonsai) < 3 || nrow(dF) < 3) {
    warning("Insufficient data: ", video)
    next
  }
  
  # ============================================================
  # 5.2 FILTER AND NORMALIZE PHOTOMETRY
  # ============================================================
  
  butter_filter <- butter(FILTER_ORDER, FILTER_CUTOFF)
  dF$dF <- as.numeric(signal::filter(butter_filter, dF$dF))
  
  max_signal <- max(dF$dF, na.rm = TRUE)
  
  if (!is.finite(max_signal) || max_signal <= 0) {
    warning("Invalid normalization maximum: ", video)
    next
  }
  
  dF$normalized_dF <- dF$dF / max_signal
  
  
  
  # ============================================================
  # 5.3 DETECT STIMULATION ONSETS
  # ============================================================
  
  episodes <- data.frame(
    episode.time = numeric(0),
    episode.number = numeric(0),
    episode.burst = numeric(0)
  )
  
  for (i in 2:(nrow(bonsai) - 1)) {
    
    if (bonsai[i, 2] > bonsai[1, 2] + PETH_WINDOW &&
        bonsai[i, 2] < bonsai[nrow(bonsai), 2] - PETH_WINDOW &&
        sum(bonsai[(i - round(PETH_WINDOW / FRAME_BEHAVIOR)):(i - 1), 3] == TRUE) == 0 &&
        bonsai[i, 3] == TRUE) {
      
      episodes[nrow(episodes) + 1, 1] <- bonsai[i, 2]
    }
  }
  
  if (nrow(episodes) == 0) {
    warning("No stimulation episodes found: ", video)
    next
  }
  
  # ============================================================
  # 5.4 IDENTIFY BURSTS AND PULSE POSITIONS
  # ============================================================
  
  episodes[1, 2] <- 1
  episodes[1, 3] <- 1
  
  if (nrow(episodes) > 1) {
    
    for (i in 2:nrow(episodes)) {
      
      if (round(episodes[i, 1]) == round(episodes[i - 1, 1] + BURST_INTERVAL)) {
        episodes[i, 2] <- as.numeric(episodes[i - 1, 2]) + 1
        episodes[i, 3] <- as.numeric(episodes[i - 1, 3])
      }
      
      if (round(episodes[i, 1]) > round(episodes[i - 1, 1] + BURST_INTERVAL)) {
        episodes[i, 2] <- 1
        episodes[i, 3] <- as.numeric(episodes[i - 1, 3]) + 1
      }
    }
  }
  
  # ============================================================
  # 5.5 EXTRACT PERI-STIMULATION SIGNALS
  # ============================================================
  
  peth_list <- vector("list", nrow(episodes))
  n_samples <- round(PETH_WINDOW / FRAME_DF)
  relative_time <- seq(-n_samples, n_samples) * FRAME_DF
  
  for (i in seq_len(nrow(episodes))) {
    
    onset <- episodes$episode.time[i]
    center <- which.min(abs(dF$time - onset))
    idx <- (center - n_samples):(center + n_samples)
    
    if (min(idx) < 1 || max(idx) > nrow(dF)) next
    
    peth_list[[i]] <- data.frame(
      nor.time = relative_time,
      normalized_dF = dF$normalized_dF[idx],
      episode.time = onset,
      episode.number = episodes$episode.number[i],
      episode.burst = episodes$episode.burst[i],
      event_id = i
    )
  }
  
  PETH <- dplyr::bind_rows(peth_list)
  
  if (nrow(PETH) == 0) {
    warning("No valid PETH windows: ", video)
    next
  }
  
  
  # ============================================================
  # 5.6 ADD EXPERIMENTAL METADATA
  # ============================================================
  
  metadata_cols <- c(
    "mouse.ID", "group", "line", "sex", "experiment",
    "day", "rec.side", "protocol", "intenisty",
    "frequency", "obs", "freq.pos", "lat.pos",
    "freq.neg", "lat.neg"
  )
  
  PETH$video <- video
  
  for (col in intersect(metadata_cols, colnames(info))) {
    PETH[[col]] <- rep(info[[col]], nrow(PETH))
  }
  
  write.csv(PETH, csv_out, row.names = FALSE)
  
  # ============================================================
  # 5.7 SUMMARIZE PETH
  # ============================================================
  
  peth_summary <- summarize_peth(PETH)
  
  # ============================================================
  # 5.8 PLOT PETH
  # ============================================================
  
 p <- ggplot(peth_summary, aes(x = time_bin, y = mean_signal)) +
  annotate("rect", xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf, fill = STIM_COLOR, alpha = 0.2) +
  geom_ribbon(aes(ymin = mean_signal - sem, ymax = mean_signal + sem),
              alpha = 0.18, na.rm = TRUE) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, color = STIM_COLOR,
             linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~ episode.number, nrow = 2) +
  scale_x_continuous(breaks = seq(-PETH_WINDOW, PETH_WINDOW, by = 1)) +
    labs(
      title = "Photometry response to optogenetic stimulation",
      subtitle = paste0("Mouse ", info$mouse.ID, " | ", info$group,
                        " | ", info$protocol, " | ", info$frequency, " Hz"),
      x = "Time from stimulation onset (s)",
      y = "Normalized fluorescence (a.u.)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.line = element_line(linewidth = 0.4),
      legend.position = "none"
    )
  
  ggsave(pdf_out, p, width = 8, height = 6)
  
  message("Completed: ", video, " | Events: ", n_distinct(PETH$event_id))
}

cat("\nPETH analysis completed.\n")
