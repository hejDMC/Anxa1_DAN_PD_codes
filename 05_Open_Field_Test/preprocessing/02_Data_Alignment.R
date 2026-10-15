# ============================================================
# VAME–DLC–fiber photometry alignment and behavioral quantification
# Author: Alessandro Contestabile
#
# Description:
#   - Aligns VAME behavioral motifs, DeepLabCut (DLC) trajectories, and fiber photometry recordings.
#   - Processes DLC trajectories and calculates locomotor speed, angular velocity, acceleration, and immobility.
#   - Applies predefined motif merging based on latent-space similarity.
#   - Calculates behavioral motif usage and summary metrics for individual recordings and mice.
#
# Requirements:
#   - details.csv containing experimental metadata.
#   - DLC tracking data exported as CSV files.
#   - VAME motif assignments stored as NumPy arrays.
#   - Fiber photometry recordings (when available).
#
# Outputs:
#   - <video>_aligned.csv: frame-wise aligned data.
#   - aligned_all.csv.gz: combined recording data.
#   - summary_per_mouse.csv: behavioral metrics and motif usage.
#
# Notes:
#   - Motif merging and preprocessing parameters are defined explicitly below to ensure reproducibility.
# ============================================================

# ------------------ CLEAN ENVIRONMENT ------------------
rm(list = ls())
invisible(gc())  

# ------------------ LIBRARIES ------------------
suppressPackageStartupMessages({
  library(dplyr)      
  library(zoo)       
  library(reticulate) 
  library(signal)     
  library(progress)   
  library(stringr)    
  library(readr)      
  library(tidyr)      
})

# ------------------ TOGGLES ------------------
OVERWRITE_ALL <- TRUE      
BUILD_BIG_ALIGNED <- TRUE  

# ------------------ PYTHON / NUMPY INIT ------------------
np <- reticulate::import("numpy", convert = TRUE)

# ------------------ DIRECTORIES  ------------------
DLC_DIR    <- "path/to/DLC_project"
VAME_DIR   <- "path/to/VAME_project"
SIGNAL_DIR <- "path/to/signals"

# Metadata (per-video) table with experimental design and flags
DETAILS_PATH <- file.path(VAME_DIR, "details.csv")

# Output root for all derived files (per-video aligned + global summaries)
OUT_ROOT <- file.path(VAME_DIR, "results")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# ------------------ PARAMETERS ------------------
# VAME model parameters
vame_k          <- 42   # number of motifs/clusters in this VAME model

# Sampling intervals (in seconds)
frame.behavior  <- 1/30 # DLC / VAME frame interval

# Threshold and scaling constants
limit           <- 0.99   # DLC body-part probability threshold for valid coordinates
px_to_cm        <- 0.0652 # DLC calibration factor: pixels → centimeters (tail base displacement)
IMMOB_THRESH    <- 0.2    # speed (cm/s) below which frames count as immobile
win             <- 30   # smoothing window (frames)

# ------------------ INPUTS ------------------
stopifnot(file.exists(DETAILS_PATH))  
details <- readr::read_csv(DETAILS_PATH, show_col_types = FALSE)

# Ensure metadata has all columns needed downstream
required_meta <- c("video","mouse.ID","group","line","sex","ctrl.vs.exp","experiment","has_signal_df")
missing_meta  <- setdiff(required_meta, names(details))
if (length(missing_meta)) stop("Missing columns in details.csv: ", paste(missing_meta, collapse = ", "))

# Normalize video IDs stored in details$video
videos <- unique(sub("\\.mp4$", "", basename(trimws(details$video)), ignore.case = TRUE))

# Detect already-aligned per-video files in OUT_ROOT
files_done <- list.files(path = OUT_ROOT, pattern = "_aligned\\.csv$", recursive = TRUE, full.names = TRUE)
processed_videos <- unique(sub("_aligned\\.csv$", "", basename(files_done)))

# Decide which videos to process in this run based on OVERWRITE_ALL toggle
videos_to_process <- if (OVERWRITE_ALL) videos else setdiff(videos, processed_videos)

message("Found ", length(videos), " videos in details; ",
        length(processed_videos), " already aligned; ",
        length(videos_to_process), " to process.",
        if (OVERWRITE_ALL) " (OVERWRITE_ALL=TRUE)" else "")

# ------------------ HELPERS ------------------
# Centered rolling mean with edge padding
smooth_vec <- function(x, k = win) {
  if (!is.numeric(x)) x <- suppressWarnings(as.numeric(x))
  if (!length(x)) return(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  y <- zoo::rollapply(
    data  = x, width = k,
    FUN   = function(xx) base::mean(xx, na.rm = TRUE),
    fill  = NA, align = "center"
  )
  # pad edges so first/last rows are not NA
  y <- zoo::na.locf(y, na.rm = FALSE)
  y <- zoo::na.locf(y, na.rm = FALSE, fromLast = TRUE)
  y
}

# Read and preprocess a single-side photometry signal (*_signal_df.csv)
read_signal <- function(file_path, butter_order = 3, butter_wc = 0.1, roll_k = win) {
  d <- readr::read_csv(file_path, show_col_types = FALSE)
  stopifnot(ncol(d) >= 2)
  names(d)[1:2] <- c("time", "dF")
  
  # Butterworth filtering (signal::filter) on raw dF
  bf <- signal::butter(butter_order, butter_wc)
  d$dF <- as.numeric(signal::filter(bf, d$dF))
  
  # Fill initial NA introduced by filtering (e.g., filter warm-up)
  d$dF <- zoo::na.locf(d$dF, na.rm = FALSE)
  
  # Apply smoothing on filtered dF
  d$dF <- smooth_vec(d$dF, k = roll_k)
  
  # Z-score the filtered and smoothed dF trace
  d$zscore <- as.numeric(scale(d$dF))
  d
}

# Nearest-neighbor alignment of behavior time stamps to signal time base
align_to_behavior <- function(behavior_time, sig_time) {
  idx <- findInterval(behavior_time, sig_time)
  idx[idx < 1] <- 1
  idx[idx > length(sig_time)] <- length(sig_time)
  idx2 <- pmin(idx + 1, length(sig_time))
  use_next <- abs(sig_time[idx2] - behavior_time) < abs(sig_time[idx] - behavior_time)
  idx[use_next] <- idx2[use_next]
  idx
}

# Find first *_signal_df.csv whose basename contains video name (case-insensitive)
find_signal_file <- function(video, signal_dir = SIGNAL_DIR) {
  all <- list.files(signal_dir, pattern = "_signal_df\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(all)) return(NA_character_)
  hit <- all[stringr::str_detect(tolower(basename(all)), stringr::fixed(tolower(video)))]
  if (!length(hit)) return(NA_character_)
  hit[1]
}

# ------------------ CORE: PROCESS ONE VIDEO ------------------
process_video <- function(video_id) {
  # ---- Select metadata row for this video ----
  meta <- details %>% dplyr::filter(.data$video == !!video_id)
  if (nrow(meta) == 0L) { 
    warning("Skipping ", video_id, " — no metadata row.")
    return(NULL)
  }
  meta <- meta[1, , drop = FALSE]
  
  # Per-video output directory and file name
  out_dir <- file.path(OUT_ROOT, video_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_csv <- file.path(out_dir, paste0(video_id, "_aligned.csv"))
  
  # Skip if already aligned and not overwriting
  if (!OVERWRITE_ALL && file.exists(out_csv)) { 
    message(video_id, " already aligned; skipping.")
    return(NULL)
  }
  
  # ---- Construct input paths for DLC and VAME labels ----
  beh_file <- file.path(DLC_DIR, paste0(video_id, ".csv"))
  vame_npy <- file.path(
    VAME_DIR, "results", video_id, paste0("VAME_", vame_k),
    paste0("hmm-", vame_k), paste0(vame_k, "_km_label_", video_id, ".npy")
  )
  
  # Existence checks for DLC and VAME
  if (!file.exists(beh_file)) { 
    warning("Skipping ", video_id, " — DLC CSV missing: ", beh_file)
    return(NULL)
  }
  if (!file.exists(vame_npy)) { 
    warning("Skipping ", video_id, " — VAME labels missing: ", vame_npy)
    return(NULL)
  }
  
  # ---- Load behavior (DLC) ----
  # DLC exports: first two rows are header/meta; data start at row 3.
  behavior <- tryCatch(
    read.csv(beh_file, skip = 2, check.names = FALSE),
    error = function(e) { 
      warning("Read error DLC for ", video_id, ": ", e$message)
      NULL
    }
  )
  if (is.null(behavior) || nrow(behavior) < 2) { 
    warning("Empty/invalid behavior for ", video_id)
    return(NULL)
  }
  
  # Standardize DLC column name scheme if columns are as expected
  expected <- c(
    "time",
    "snout_x","snout_y","snout_prob",
    "leftForelimb_x","leftForelimb_y","leftForelimb_prob",
    "rightForelimb_x","rightForelimb_y","rightForelimb_prob",
    "leftPosteriorlimb_x","leftPosteriorlimb_y","leftPosteriorlimb_prob",
    "rightPosteriorlimb_x","rightPosteriorlimb_y","rightPosteriorlimb_prob",
    "tailBase_x","tailBase_y","tailBase_prob"
  )
  if (ncol(behavior) >= length(expected)) {
    colnames(behavior)[seq_along(expected)] <- expected
  } else {
    warning("Non-standard DLC columns for ", video_id, " — continuing with current names.")
  }
  
  # Convert frame index to seconds using behavior frame interval
  behavior$time <- behavior$time * frame.behavior
  
  # ---- Load VAME motifs ----
  vm <- tryCatch(
    np$load(vame_npy),
    error = function(e) { 
      warning("VAME load error for ", video_id, ": ", e$message)
      NULL
    }
  )
  if (is.null(vm) || length(vm) < 2) { 
    warning("Invalid VAME labels for ", video_id)
    return(NULL)
  }
  vame_motifs <- data.frame(motif = as.integer(vm))
  
  # ---- Match behavior and VAME lengths ----
  n_use <- min(nrow(behavior), nrow(vame_motifs))
  behavior    <- behavior   [seq_len(n_use), , drop = FALSE]
  vame_motifs <- vame_motifs[seq_len(n_use), , drop = FALSE]
  behavior$motif <- vame_motifs$motif
  
  # ---- Motif merging ----
  # Predefined motif merges based on inspection of latent-space similarity (see VAME dendrogram analysis).
  behavior$motif[behavior$motif == 11] <- 28
  behavior$motif[behavior$motif == 2] <- 8
  behavior$motif[behavior$motif == 15] <- 38
  behavior$motif[behavior$motif == 29] <- 34
  
  # ---- Clean / interpolate DLC points ----
  # Identify body-part roots (snout, limbs, tailBase, etc.) from *_x/_y/_prob
  point_roots <- unique(gsub("_(x|y|prob)$", "", names(behavior)))
  point_roots <- setdiff(point_roots, c("time","motif"))
  
  for (p in point_roots) {
    px <- paste0(p, "_x")
    py <- paste0(p, "_y")
    pp <- paste0(p, "_prob")
    if (!all(c(px,py,pp) %in% names(behavior))) next
    
    # Drop low-confidence coordinates (prob < limit) to NA
    behavior[[px]] <- ifelse(behavior[[pp]] < limit, NA, behavior[[px]])
    behavior[[py]] <- ifelse(behavior[[pp]] < limit, NA, behavior[[py]])
    
    # Interpolate missing coordinates over time (linear, with rule=2 for edges)
    behavior[[px]] <- zoo::na.approx(behavior[[px]], na.rm = FALSE, rule = 2)
    behavior[[py]] <- zoo::na.approx(behavior[[py]], na.rm = FALSE, rule = 2)
    
    # Smooth x/y trajectories with rolling mean (same window as used elsewhere)
    behavior[[px]] <- smooth_vec(behavior[[px]], k = win)
    behavior[[py]] <- smooth_vec(behavior[[py]], k = win)
  }
  
  # ---- Angular and linear speed ----
  # Angular speed: orientation of snout–tailBase axis over time
  if (all(c("snout_x","snout_y","tailBase_x","tailBase_y") %in% names(behavior))) {
    # vectors from tailBase to snout at time t and t+1
    ax <- behavior$snout_x[-nrow(behavior)] - behavior$tailBase_x[-nrow(behavior)]
    ay <- behavior$snout_y[-nrow(behavior)] - behavior$tailBase_y[-nrow(behavior)]
    bx <- behavior$snout_x[-1]              - behavior$tailBase_x[-1]
    by <- behavior$snout_y[-1]              - behavior$tailBase_y[-1]
    
    # cross/dot product to compute signed angle between successive vectors
    cross <- ax * by - ay * bx
    dot   <- ax * bx + ay * by
    ang   <- atan2(cross, dot)                         
    ang[!is.finite(ang)] <- NA_real_
    
    # Angular distance per frame and angular speed (|Δθ| / Δt)
    behavior$angular.distance <- c(NA_real_, abs(ang))
    behavior$angular.speed    <- behavior$angular.distance / frame.behavior
  } else {
    behavior$angular.distance <- NA_real_
    behavior$angular.speed    <- NA_real_
  }
  
  # Snout speed (cm/s)
  if (all(c("snout_x","snout_y") %in% names(behavior))) {
    dx_snout <- diff(behavior$snout_x)
    dy_snout <- diff(behavior$snout_y)
    behavior$distance.moved_snout <- c(NA_real_, sqrt(dx_snout^2 + dy_snout^2) * px_to_cm)
    behavior$speed_snout          <- behavior$distance.moved_snout / frame.behavior
  } else {
    behavior$distance.moved_snout <- NA_real_
    behavior$speed_snout          <- NA_real_
  }
  
  # Tail-base speed (cm/s)
  if (all(c("tailBase_x","tailBase_y") %in% names(behavior))) {
    dx_tail <- diff(behavior$tailBase_x)
    dy_tail <- diff(behavior$tailBase_y)
    behavior$distance.moved <- c(NA_real_, sqrt(dx_tail^2 + dy_tail^2) * px_to_cm)
    behavior$speed          <- behavior$distance.moved / frame.behavior
    
    # Linear acceleration (cm/s^2) from tail-base speed
    behavior$acceleration   <- c(NA_real_, diff(behavior$speed)) / frame.behavior
  } else {
    behavior$distance.moved <- NA_real_
    behavior$speed          <- NA_real_
    behavior$acceleration   <- NA_real_
  }
  
  # ---- Photometry (single side: *_signal_df.csv) ----
  # has_signal_df column determines whether photometry is expected for this video
  has_sig <- as.logical(toupper(as.character(meta$has_signal_df)) %in% c("TRUE"))
  behavior$signal <- NA_real_
  
  if (has_sig) {
    sig_file <- find_signal_file(video_id, SIGNAL_DIR)
    if (!is.na(sig_file) && file.exists(sig_file)) {
      sig <- tryCatch(read_signal(sig_file), error = function(e) NULL)
      if (!is.null(sig) && nrow(sig) > 1) {
        # Align photometry z-score to behavior time base via nearest neighbor
        idx <- align_to_behavior(behavior$time, sig$time)
        behavior$signal <- sig$zscore[idx]
      } else {
        message("Signal unreadable/empty for ", video_id, " — continuing without signal.")
      }
    } else {
      message("No *_signal_df.csv matched for ", video_id, " — continuing without signal.")
    }
  }
  
  # ---- Attach metadata ----
  behavior$mouse.ID     <- as.character(meta$mouse.ID)
  behavior$group        <- as.character(meta$group)
  behavior$line         <- as.character(meta$line)
  behavior$sex          <- as.character(meta$sex)
  behavior$ctrl.vs.exp  <- as.character(meta$ctrl.vs.exp)
  behavior$experiment   <- as.character(meta$experiment)
  behavior$video        <- video_id
  
  # ---- SAVE per-video aligned ----
  readr::write_csv(behavior, out_csv, na = "")
  
  # ---- Per-video summary (returned as one row per video) ----
  mean_speed    <- mean(behavior$speed, na.rm = TRUE)
  mean_ang_spd  <- mean(behavior$angular.speed, na.rm = TRUE)
  
  # immobility: require BOTH linear speeds below IMMOB_THRESH when available
  valid_tail  <- !is.na(behavior$speed)
  valid_snout <- !is.na(behavior$speed_snout)
  mean_accel <- mean(behavior$acceleration, na.rm = TRUE)
  
  if (any(valid_snout)) {
    # both speeds available → define immobility on frames where both are valid
    valid_both    <- valid_tail & valid_snout
    immobile_rows <- valid_both &
      (behavior$speed < IMMOB_THRESH) &
      (behavior$speed_snout < IMMOB_THRESH)
    
    immobile_time <- sum(immobile_rows) * frame.behavior
    total_time_s  <- sum(valid_both) * frame.behavior
  } else {
    immobile_rows <- valid_tail & (behavior$speed < IMMOB_THRESH)
    immobile_time <- sum(immobile_rows) * frame.behavior
    total_time_s  <- sum(valid_tail) * frame.behavior
  }
  
  immobile_pct  <- if (total_time_s > 0) 100 * immobile_time / total_time_s else NA_real_
  
  # Motif usage: total time and percentage per motif (within this video)
  valid_motif <- !is.na(behavior$motif)
  motif_tbl   <- table(behavior$motif[valid_motif])
  motif_sec   <- as.numeric(motif_tbl) * frame.behavior
  motif_pct   <- if (sum(valid_motif) > 0) {
    100 * as.numeric(motif_tbl) / sum(valid_motif)
  } else {
    rep(NA_real_, length(motif_tbl))
  }
  
  motif_df <- tibble::tibble(
    motif     = as.integer(names(motif_tbl)),
    motif_sec = motif_sec,
    motif_pct = motif_pct
  )
  
  # Wide representation: one column per motif (seconds and percent)
  motif_sec_w <- tidyr::pivot_wider(
    motif_df, names_from = motif, values_from = motif_sec, names_prefix = "motif_sec_"
  )
  motif_pct_w <- tidyr::pivot_wider(
    motif_df, names_from = motif, values_from = motif_pct, names_prefix = "motif_pct_"
  )
  
  video_summary <- tibble::tibble(
    video           = video_id,
    mouse.ID        = behavior$mouse.ID[1],
    group           = behavior$group[1],
    line            = behavior$line[1],
    sex             = behavior$sex[1],
    ctrl.vs.exp     = behavior$ctrl.vs.exp[1],
    experiment      = behavior$experiment[1],
    mean_speed      = mean_speed,
    mean_ang_spd    = mean_ang_spd,
    mean_accel      = mean_accel,
    immobile_time_s = immobile_time,
    immobile_pct    = immobile_pct,
    total_time_s    = total_time_s
  ) %>% 
    dplyr::bind_cols(motif_sec_w) %>% 
    dplyr::bind_cols(motif_pct_w)
  
  return(video_summary)
}

# ------------------ RUN ------------------
summaries <- list()

# Progress bar reports video index and ID (token :video)
pb <- progress::progress_bar$new(
  format = "[:bar] :current/:total :percent :eta | :video",
  total = length(videos_to_process), clear = FALSE, width = 80
)

for (vid in videos_to_process) {
  res <- tryCatch(
    process_video(vid),
    error = function(e) {
      warning("Failed processing ", vid, ": ", e$message)
      NULL
    }
  )
  if (!is.null(res)) summaries[[length(summaries) + 1L]] <- res
  pb$tick(tokens = list(video = vid))
}

message("Processed videos (this run): ", length(summaries))

# ------------------ BIG ALIGNED FILE ------------------
big_aligned_path <- file.path(OUT_ROOT, "aligned_all.csv.gz")

if (BUILD_BIG_ALIGNED) {
  # Columns that should always be character (metadata)
  meta_char_cols <- c("video","mouse.ID","group","line","sex","ctrl.vs.exp","experiment")
  
  # Readr column specification to enforce types for critical metadata columns
  col_spec <- readr::cols(
    .default     = readr::col_guess(),
    video        = readr::col_character(),
    mouse.ID     = readr::col_character(),
    group        = readr::col_character(),
    line         = readr::col_character(),
    sex          = readr::col_character(),
    ctrl.vs.exp  = readr::col_character(),
    experiment   = readr::col_character()
  )
  
  # Safe reader that enforces metadata types and numeric casts where appropriate
  safe_read_aligned <- function(path) {
    df <- tryCatch(
      readr::read_csv(path, col_types = col_spec, show_col_types = FALSE),
      error = function(e) NULL
    )
    if (is.null(df)) return(NULL)
    
    for (cc in meta_char_cols) {
      if (!cc %in% names(df)) df[[cc]] <- NA_character_
      df[[cc]] <- as.character(df[[cc]])
    }
    
    num_patterns <- c("(_x$|_y$)", "^speed$", "^angular\\.speed$", "^distance\\.moved$", "^signal$", "^time$", "^acceleration$")
    num_cols <- grep(paste(num_patterns, collapse="|"), names(df), value = TRUE)
    for (nc in num_cols) df[[nc]] <- suppressWarnings(as.numeric(df[[nc]]))
    
    if ("motif" %in% names(df)) df$motif <- suppressWarnings(as.integer(df$motif))
    df
  }
  
  # Gather all per-video aligned CSVs
  aligned_csvs <- list.files(OUT_ROOT, pattern = "_aligned\\.csv$", recursive = TRUE, full.names = TRUE)
  aligned_list <- lapply(aligned_csvs, safe_read_aligned)
  aligned_list <- aligned_list[!vapply(aligned_list, is.null, logical(1))]
  
  if (length(aligned_list)) {
    big_aligned <- dplyr::bind_rows(aligned_list)
    
    # Edge padding: remove leading/trailing gaps in key trajectories and signals
    fill_cols <- grep("(_x$|_y$|^speed$|^angular\\.speed$|^distance\\.moved$|^signal$|^acceleration$)",
                      names(big_aligned), value = TRUE)
    if (length(fill_cols)) {
      big_aligned[fill_cols] <- lapply(big_aligned[fill_cols], function(col) {
        col <- zoo::na.locf(col, na.rm = FALSE)
        col <- zoo::na.locf(col, na.rm = FALSE, fromLast = TRUE)
        col
      })
    }
    
    # Write gzipped CSV for downstream use (e.g., PETHs, motif analyses, etc.)
    con <- gzfile(big_aligned_path, open = "wb")
    on.exit(close(con), add = TRUE)
    readr::write_csv(big_aligned, con)
    message("Wrote concatenated aligned file: ", big_aligned_path,
            " (rows: ", nrow(big_aligned), ")")
  } else {
    message("No aligned CSVs found to concatenate.")
  }
}

# ------------------ BUILD PER-MOUSE SUMMARY ------------------
make_per_mouse <- function(df_all) {
  keys <- c("mouse.ID","group","line","sex","ctrl.vs.exp","experiment")
  
  per_mouse_core <- df_all %>%
    dplyr::group_by(dplyr::across(all_of(keys))) %>%
    dplyr::summarise(
      n_videos        = dplyr::n_distinct(.data$video),
      total_time_s    = sum(!is.na(.data$speed)) * frame.behavior,
      mean_speed      = base::mean(.data$speed, na.rm = TRUE),
      mean_ang_speed  = base::mean(.data$angular.speed, na.rm = TRUE),
      mean_accel      = base::mean(.data$acceleration, na.rm = TRUE),
      immobile_time_s = sum(!is.na(.data$speed) & .data$speed < IMMOB_THRESH) * frame.behavior,
      .groups = "drop"
    )
  
  motif_counts <- df_all %>%
    dplyr::filter(!is.na(.data$motif)) %>%
    dplyr::group_by(dplyr::across(all_of(keys)), .data$motif) %>%
    dplyr::summarise(n_frames = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(motif_sec = .data$n_frames * frame.behavior)
  
  motif_sec_wide <- motif_counts %>%
    dplyr::select(-.data$n_frames) %>%
    tidyr::pivot_wider(
      names_from = .data$motif,
      values_from = .data$motif_sec,
      values_fill = 0,
      names_prefix = "motif_sec_"
    )
  
  motif_sec_cols <- grep("^motif_sec_", names(motif_sec_wide), value = TRUE)
  motif_pct_wide <- motif_sec_wide
  if (length(motif_sec_cols)) {
    motif_pct_wide <- motif_sec_wide %>%
      dplyr::rowwise() %>%
      dplyr::mutate(total_motif_sec = sum(c_across(dplyr::all_of(motif_sec_cols)), na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(motif_sec_cols),
          ~ dplyr::if_else(total_motif_sec > 0, 100 * .x / total_motif_sec, NA_real_),
          .names = "{gsub('motif_sec_', 'motif_pct_', .col)}"
        )
      ) %>%
      dplyr::select(dplyr::all_of(keys), dplyr::starts_with("motif_pct_"))
  } else {
  
    motif_pct_wide <- dplyr::select(per_mouse_core, dplyr::all_of(keys))
  }
  
  per_mouse_core %>%
    dplyr::left_join(motif_sec_wide, by = keys) %>%
    dplyr::left_join(motif_pct_wide, by = keys) %>%
    dplyr::arrange(.data$group, .data$mouse.ID)
}

per_mouse_full <- NULL

# Prefer reading from aligned_all.csv.gz (if it exists and BUILD_BIG_ALIGNED is TRUE)
if (BUILD_BIG_ALIGNED && file.exists(big_aligned_path)) {
  df_all <- readr::read_csv(big_aligned_path, show_col_types = FALSE)
  
  # Ensure expected types (esp. motif and metadata) even after external edits
  if ("motif" %in% names(df_all)) df_all$motif <- suppressWarnings(as.integer(df_all$motif))
  meta_char_cols <- c("video","mouse.ID","group","line","sex","ctrl.vs.exp","experiment")
  for (cc in meta_char_cols) if (cc %in% names(df_all)) df_all[[cc]] <- as.character(df_all[[cc]])
  
  per_mouse_full <- make_per_mouse(df_all)
  
} else if (length(summaries)) {
  aligned_csvs <- list.files(OUT_ROOT, pattern = "_aligned\\.csv$", recursive = TRUE, full.names = TRUE)
  aligned_list <- lapply(aligned_csvs, function(f) {
    tryCatch(readr::read_csv(f, show_col_types = FALSE), error = function(e) NULL)
  })
  aligned_list <- aligned_list[!vapply(aligned_list, is.null, logical(1))]
  if (length(aligned_list)) {
    df_all <- dplyr::bind_rows(aligned_list)
    if ("motif" %in% names(df_all)) df_all$motif <- suppressWarnings(as.integer(df_all$motif))
    per_mouse_full <- make_per_mouse(df_all)
  }
}

# Write per-mouse summary if available
if (!is.null(per_mouse_full)) {
  out_summary <- file.path(OUT_ROOT, "summary_per_mouse.csv")
  readr::write_csv(per_mouse_full, out_summary, na = "")
  message("Wrote per-mouse summary: ", out_summary,
          "  (rows=", nrow(per_mouse_full),
          "; mice=", dplyr::n_distinct(per_mouse_full$mouse.ID), ")")
} else {
  message("No summaries to write.")
}

# ------------------ RUN RECAP ------------------
if (length(summaries)) {
  u_mice <- unique(dplyr::bind_rows(summaries)$mouse.ID)
  message("Mice covered this run: ", length(u_mice), " [",
          paste(utils::head(u_mice, 8), collapse = ", "),
          if (length(u_mice) > 8) ", …" else "", "]")
}