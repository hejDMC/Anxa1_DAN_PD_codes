# ============================================================
# Figure 5D–H — Open-field behavioral motifs and photometry
# Author: Alessandro Contestabile
#
# Description:
#   - Quantifies behavioral motif duration and locomotor speed.
#   - Characterizes motif distributions using ridgeline plots.
#   - Generates motif-aligned photometry heatmaps.
#   - Calculates post-onset photometry responses (AUC, 0–1 s).
#   - Examines relationships between motif speed, dopamine signals, and Anxa1+ neuron activity.
#   - Assembles the Figure 5D–F panels.
#
# Requirements:
#   - Frame-wise aligned VAME, DLC, and photometry data.
#   - Motif-aligned PETH summaries from preprocessing.
#
# Outputs:
#   - figure_5D-H.pdf
#   - Motif duration and speed summaries.
#   - Motif-aligned photometry heatmaps.
#   - AUC summaries and supplementary figures.
# ============================================================

# ------------------ CLEAN ENV ------------------
rm(list = ls()); invisible(gc())

# ------------------ LIBRARIES ------------------
suppressPackageStartupMessages({
  library(tidyverse)   
  library(ggridges)   
  library(ggpubr)     
  library(patchwork)   
  library(zoo)   
})

# ------------------ GLOBAL SETTINGS ------------------
set.seed(1234)

# Root containing *_aligned.csv and *_peth_agg.csv
OUT_ROOT <- "path/to/VAME_project/results"

# Analysis directory for all outputs
ANALYSIS_DIR <- file.path(OUT_ROOT, "analysis")
dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)

# Save session info for reproducibility
writeLines(c(
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "------------- SESSION INFO -------------",
  capture.output(sessionInfo())
), con = file.path(ANALYSIS_DIR, "session_info.txt"))

# Behavior sampling interval (s) — used for motif episode duration
frame.behavior <- 1/30

# Motif inclusion criteria for all graphs
MIN_MICE_PER_MOTIF <- 40

# Videos to exclude entirely from all analyses (if any)
# e.g. VIDEO_TO_REMOVE <- c("oft_2025-11-05T15_58_07_Left", "...")
VIDEO_TO_REMOVE <- c("oft_2023-11-28T14_02_24", "oft_2023-11-28T14_28_05", "oft_2023-11-28T13_36_09", "oft_2024-06-05T14_53_19", "tmaze_2023-05-15T14_13_39")

# PETH binning settings
BIN_WIDTH <- 0.05                   # 50 ms bins
XLIM      <- c(-1.5, 3)             # time window (s) for PETH plots
PREFERRED_SIGNALS  <- c()   # if empty, use all signals
INVERT_FACET_ORDER <- FALSE

# Manual motif exclusions for PETH (numeric IDs) if needed
EXCLUDE_MOTIFS <- integer(0)        # e.g. c(21, 29, 10)

# For heatmap filenames
suffix <- sprintf("%.0fms", BIN_WIDTH * 1000)

# ------------------ HELPER FUNCTIONS ------------------
aligned_col_spec <- readr::cols(
  .default     = readr::col_guess(),
  time         = readr::col_double(),
  video        = readr::col_character(),
  mouse.ID     = readr::col_character(),
  group        = readr::col_character(),
  line         = readr::col_character(),
  sex          = readr::col_character(),
  ctrl.vs.exp  = readr::col_character(),
  experiment   = readr::col_character(),
  motif        = readr::col_integer(),
  speed        = readr::col_double()
)

safe_read_aligned <- function(path) {
  tryCatch(
    readr::read_csv(path, col_types = aligned_col_spec, show_col_types = FALSE),
    error = function(e) NULL
  )
}

num_coerce    <- function(x) suppressWarnings(as.numeric(x))
safe_filename <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
bin_time      <- function(x, width = BIN_WIDTH) floor(x / width) * width

read_csv_safely_base <- function(f) {
  tryCatch(
    read.csv(f, check.names = FALSE),
    error = function(e) { message("!! Failed to read: ", f); NULL }
  )
}

# ------------------ 1) OVERALL MEAN SPEED PER MOTIF ------------------
aligned_csvs <- list.files(
  OUT_ROOT,
  pattern   = "_aligned\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(aligned_csvs)) {
  stop("No '*_aligned.csv' files found under OUT_ROOT: ", OUT_ROOT)
}

aligned_list <- lapply(aligned_csvs, safe_read_aligned)
aligned_list <- aligned_list[!vapply(aligned_list, is.null, logical(1))]
if (!length(aligned_list)) stop("All aligned files failed to read or are empty.")

dat <- dplyr::bind_rows(aligned_list)

req_cols <- c("motif","speed","video","mouse.ID","experiment","time")
missing_cols <- setdiff(req_cols, names(dat))
if (length(missing_cols)) {
  stop("Missing required columns in aligned data: ", paste(missing_cols, collapse = ", "))
}

# Remove selected videos globally
if (length(VIDEO_TO_REMOVE) && "video" %in% names(dat)) {
  dat <- dat %>% dplyr::filter(!(video %in% VIDEO_TO_REMOVE))
}

# Experiments contributing to motif speed ordering
ALLOWED_EXPERIMENTS <- c(
  "oft_6OHDA",
  "oft_mitopark_10-11weeks",
  "oft_mitopark_15-18weeks",
  "oft_Tetx",
  "oft_fp"
)

dat_filtered <- dat %>%
  dplyr::filter(experiment %in% ALLOWED_EXPERIMENTS)

# Pooled mean speed per motif (across all videos & mice)
motif_mean_overall <- dat_filtered %>%
  dplyr::filter(!is.na(motif), !is.na(speed)) %>%
  dplyr::group_by(motif) %>%
  dplyr::summarise(
    mean_speed_motif = mean(speed, na.rm = TRUE),
    sd_motif         = sd(speed,   na.rm = TRUE),
    n_frames         = dplyr::n(),
    n_videos         = dplyr::n_distinct(video),
    n_mice           = dplyr::n_distinct(mouse.ID),
    .groups          = "drop"
  ) %>%
  dplyr::arrange(motif)

out_overall <- file.path(ANALYSIS_DIR, "mean_speed_per_motif_overall.csv")
readr::write_csv(motif_mean_overall, out_overall)
message("Saved: ", out_overall, " (rows: ", nrow(motif_mean_overall), ")")

# Motifs with sufficient representation across mice
valid_motifs <- motif_mean_overall %>%
  dplyr::filter(n_mice >= MIN_MICE_PER_MOTIF) %>%
  dplyr::pull(motif)

message("Motifs retained for graphs (n_mice ≥ ", MIN_MICE_PER_MOTIF, "): ",
        length(valid_motifs))

# ------------------ 2) MOTIF EPISODES: LENGTH and EPISODE MEAN SPEED ------------------
episodes <- dat_filtered %>%
  dplyr::filter(!is.na(motif), !is.na(time), motif %in% valid_motifs) %>%
  dplyr::arrange(video, time) %>%
  dplyr::group_by(video, mouse.ID, experiment, group) %>%
  dplyr::group_modify(function(df, key) {
    m <- df$motif
    if (!length(m)) return(tibble())
    
    r   <- rle(m)
    len <- r$lengths
    val <- r$values
    end_idx   <- cumsum(len)
    start_idx <- c(1L, head(end_idx, -1) + 1L)
    
    start_time <- df$time[start_idx]
    end_time   <- df$time[end_idx]
    motif_len  <- pmax(0, (end_time - start_time) + frame.behavior)
    
    mean_speed_ep <- purrr::map2_dbl(
      start_idx, end_idx,
      ~ mean(df$speed[.x:.y], na.rm = TRUE)
    )
    
    tibble(
      motif         = val,
      start_idx     = start_idx,
      end_idx       = end_idx,
      start_time    = start_time,
      end_time      = end_time,
      motif_length  = motif_len,
      mean_speed_ep = mean_speed_ep
    )
  }) %>%
  dplyr::ungroup()

out_episodes <- file.path(ANALYSIS_DIR, "motif_episodes.csv")
readr::write_csv(episodes, out_episodes)
message("Saved: ", out_episodes, " (rows: ", nrow(episodes), ")")

# ------------------ 3) RIDGELINE PLOTS: MOTIF LENGTH & EPISODE SPEED ------------------
# Attach global mean_speed_motif to each episode
episodes <- episodes %>%
  dplyr::left_join(
    motif_mean_overall %>%
      dplyr::select(motif, mean_speed_motif),
    by = "motif"
  )

# Order motifs by mean_speed_motif (slow → fast; fast on top in plots)
motif_levels <- episodes %>%
  dplyr::group_by(motif) %>%
  dplyr::summarise(
    mean_speed_motif = dplyr::first(stats::na.omit(mean_speed_motif)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mean_speed_motif) %>%
  dplyr::pull(motif) %>%
  as.character()

episodes <- episodes %>%
  dplyr::mutate(
    motif_f     = factor(as.character(motif), levels = motif_levels),
    motif_index = as.numeric(motif_f)  # numeric y-position per motif
  )

# Median motif length per motif (for vertical ticks)
med_len <- episodes %>%
  dplyr::group_by(motif_f, motif_index, mean_speed_motif) %>%
  dplyr::summarise(
    median_len = stats::median(motif_length, na.rm = TRUE),
    .groups    = "drop"
  )

# Median episode speed per motif
med_speed <- episodes %>%
  dplyr::group_by(motif_f, motif_index, mean_speed_motif) %>%
  dplyr::summarise(
    median_speed = stats::median(mean_speed_ep, na.rm = TRUE),
    .groups      = "drop"
  )

# ---- Motif length ridgelines ----
p_len <- ggplot(
  episodes,
  aes(
    x     = motif_length,
    y     = motif_f,
    group = motif_f,
    fill  = mean_speed_motif
  )
) +
  ggridges::geom_density_ridges(scale = 5, alpha = 0.9, color = "black") +
  geom_segment(
    data = med_len,
    aes(
      x    = median_len,
      xend = median_len,
      y    = motif_index + 1,
      yend = motif_index
    ),
    inherit.aes = FALSE,
    linewidth   = 1,
    color       = "black"
  ) +
  coord_cartesian(xlim = c(0, 4)) +
  scale_fill_distiller(
    palette   = "BuPu",
    direction = 1,
    name      = "Mean speed\nper motif"
  ) +
  xlab("Motif length (s)") +
  ylab("Motif") +
  theme_bw() +
  theme(
    legend.position   = "none",
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.background  = element_blank()
  )

# ---- Episode mean-speed ridgelines ----
p_speed <- ggplot(
  episodes,
  aes(
    x     = mean_speed_ep,
    y     = motif_f,
    group = motif_f,
    fill  = mean_speed_motif
  )
) +
  ggridges::geom_density_ridges(scale = 10, alpha = 0.9, color = "black") +
  geom_segment(
    data = med_speed,
    aes(
      x    = median_speed,
      xend = median_speed,
      y    = motif_index + 1,
      yend = motif_index
    ),
    inherit.aes = FALSE,
    linewidth   = 1,
    color       = "black"
  ) +
  coord_cartesian(xlim = c(0, 20)) +
  scale_fill_distiller(
    palette   = "BuPu",
    direction = 1,
    name      = "Mean speed\nper motif"
  ) +
  xlab("Mean speed per motif episode (cm/s)") +
  ylab("Motif") +
  theme_bw() +
  theme(
    legend.position   = "none",
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.background  = element_blank()
  )

# ------------------ 4) MOTIF-ALIGNED PETH HEATMAPS ORDERED BY MEAN SPEED ------------------

# ------------------ LOAD ALL *_peth_agg.csv ------------------
peth_files <- list.files(
  path      = OUT_ROOT,
  pattern   = "_peth_agg.csv",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(peth_files)) {
  stop("No '*_peth_agg.csv' files found under OUT_ROOT: ", OUT_ROOT)
}

combined_PETH <- purrr::map_dfr(
  peth_files,
  \(f) {
    df <- read_csv_safely_base(f)
    if (is.null(df)) return(NULL)
    dplyr::mutate(df, source_file = f)
  }
)

req_peth_cols <- c("event_var", "signal", "n_events", "mean_value", "nor.time")
missing_peth  <- setdiff(req_peth_cols, names(combined_PETH))
if (length(missing_peth)) {
  stop("combined_PETH missing cols: ", paste(missing_peth, collapse = ", "))
}

# Drop selected videos from PETH as well
if (length(VIDEO_TO_REMOVE) && "video" %in% names(combined_PETH)) {
  combined_PETH <- combined_PETH %>%
    dplyr::filter(!(video %in% VIDEO_TO_REMOVE))
}

if (!"group" %in% names(combined_PETH)) {
  stop("Column 'group' not found in *_peth_agg.csv; cannot split plots by group.")
}

# ------------------ HOW MANY VIDEOS PER GROUP (PETH) ------------------
group_summary <- combined_PETH %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    n_videos = dplyr::n_distinct(video),
    .groups  = "drop"
  ) %>%
  dplyr::arrange(group)

# Named vector for title annotation
n_videos_title <- stats::setNames(group_summary$n_videos, group_summary$group)

# Save table
out_peth_group <- file.path(ANALYSIS_DIR, "peth_videos_per_group.csv")
readr::write_csv(group_summary, out_peth_group)
message("Saved PETH contributors per group (videos): ", out_peth_group)

# ------------------ KEEP ONLY MOTIF EVENTS ------------------
peth_motif <- combined_PETH %>%
  dplyr::filter(grepl("^\\s*motif[_ ]*\\d+\\s*$", event_var)) %>%
  dplyr::mutate(
    motif_num  = suppressWarnings(as.integer(gsub(".*?(\\d+)\\s*$", "\\1", event_var))),
    n_events   = num_coerce(n_events),
    mean_value = num_coerce(mean_value),
    sem_value  = num_coerce(sem_value),
    nor.time   = num_coerce(nor.time),
    group      = as.character(group),
    signal     = as.character(signal)
  ) %>%
  dplyr::filter(
    !is.na(motif_num),
    !is.na(nor.time),
    motif_num %in% valid_motifs,
    !(motif_num %in% EXCLUDE_MOTIFS)
  )

signals_all     <- sort(unique(peth_motif$signal))
signals_in_plot <- PREFERRED_SIGNALS[PREFERRED_SIGNALS %in% signals_all]
if (!length(signals_in_plot)) signals_in_plot <- signals_all
if (INVERT_FACET_ORDER) signals_in_plot <- rev(signals_in_plot)

groups_all <- sort(unique(peth_motif$group))

# ------------------ KEEP ONLY SPEED PETH ------------------
peth_speed <- combined_PETH %>%
  dplyr::filter(grepl("^\\s*motif[_ ]*\\d+\\s*$", event_var), signal == "speed") %>%
  dplyr::mutate(
    motif_num  = suppressWarnings(as.integer(gsub(".*?(\\d+)\\s*$", "\\1", event_var))),
    n_events   = num_coerce(n_events),
    mean_value = num_coerce(mean_value),
    sem_value  = num_coerce(sem_value),
    nor.time   = num_coerce(nor.time),
    group      = as.character(group),
    signal     = as.character(signal)
  ) %>%
  dplyr::filter(
    !is.na(motif_num),
    !is.na(nor.time),
    motif_num %in% valid_motifs,
    !(motif_num %in% EXCLUDE_MOTIFS)
  )

# ------------------ SPEED-BASED ORDER (DIRECTLY FROM motif_mean_overall) ------------------
ord_tbl <- motif_mean_overall %>%
  dplyr::filter(motif %in% valid_motifs, !(motif %in% EXCLUDE_MOTIFS)) %>%
  dplyr::transmute(
    motif_num        = as.integer(motif),
    mean_speed_motif = num_coerce(mean_speed_motif)
  ) %>%
  dplyr::filter(!is.na(motif_num))

# ------------------ BIN & MERGE MEAN SPEED INTO PETH ------------------
agg_motif_binned <- peth_motif %>%
  dplyr::mutate(nor_bin = bin_time(nor.time, BIN_WIDTH)) %>%
  dplyr::group_by(group, motif_num, signal, nor_bin) %>%
  dplyr::summarise(
    # total n across sub-bins
    n_events_total = sum(n_events, na.rm = TRUE),
    
    # pooled mean
    mean_value_pooled = ifelse(
      n_events_total > 0,
      sum(mean_value * n_events, na.rm = TRUE) / n_events_total,
      NA_real_
    ),
    
    # recover per-row variance from SEM: var_i = SEM_i^2 * n_i
    pooled_var = dplyr::if_else(
      n_events_total > 1,
      {
        row_var <- (sem_value ^ 2) * n_events
        mu      <- sum(mean_value * n_events, na.rm = TRUE) / n_events_total
        
        # classic pooled-variance formula
        num <- sum(
          (n_events - 1) * row_var +
            n_events * (mean_value - mu) ^ 2,
          na.rm = TRUE
        )
        
        num / (n_events_total - 1)
      },
      NA_real_
    ),
    
    sd_value  = sqrt(pooled_var),
    sem_value = sd_value / sqrt(n_events_total),
    
    # final mean we keep
    mean_value = mean_value_pooled,
    
    .groups = "drop"
  ) %>%
  dplyr::filter(nor_bin >= XLIM[1], nor_bin <= XLIM[2])


speed_info <- ord_tbl %>%
  dplyr::select(motif_num, mean_speed_motif) %>%
  dplyr::distinct()

agg_motif_binned <- agg_motif_binned %>%
  dplyr::left_join(speed_info, by = "motif_num") %>%
  dplyr::mutate(sort_speed = ifelse(is.na(mean_speed_motif), -Inf, mean_speed_motif))

# ------------------ PETH HEATMAPS PER GROUP ------------------
peth_plots_full <- list()

for (g in groups_all) {
  df_g <- agg_motif_binned %>%
    dplyr::filter(group == g, signal == "signal") %>%
    dplyr::arrange(dplyr::desc(sort_speed), motif_num)
  
  if (!nrow(df_g)) next
  
  motif_levels_g <- unique(df_g$motif_num)
  df_g <- df_g %>%
    dplyr::mutate(
      motif_num = factor(motif_num, levels = motif_levels_g),
      signal    = factor(signal, levels = signals_in_plot)
    )
  
  # group-level n_videos for title
  n_text <- ""
  if (length(n_videos_title) && g %in% names(n_videos_title)) {
    n_text <- sprintf(" (n = %d videos)", n_videos_title[[g]])
  }
  
  # Full plot (with axes and legend)
  p_g_base <- ggplot(df_g, aes(x = nor_bin, y = motif_num, fill = mean_value)) +
    geom_tile() +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_fill_distiller(
      type    = "div",
      palette = "PuOr",
      limits  = c(-0.3, 0.3),
      oob     = scales::squish
    ) +
    coord_cartesian(xlim = XLIM, expand = FALSE) +
    scale_y_discrete(limits = rev(levels(df_g$motif_num))) +   # fastest on top
    labs(
      title = paste0(g, n_text),
      x     = "Time (s)",
      y     = "Motif",
      fill  = "Mean z"
    ) +
    theme_bw() +
    theme(
      panel.grid      = element_blank(),
      legend.position = "bottom"
    )
  
  peth_plots_full[[g]] <- p_g_base
  
  # Minimal heatmap (no axes; title only) for PNG export
  p_g_heat <- ggplot(df_g, aes(x = nor_bin, y = motif_num, fill = mean_value)) +
    geom_tile() +
    scale_fill_distiller(
      type    = "div",
      palette = "PuOr",
      limits  = c(-0.3, 0.3),
      oob     = scales::squish
    ) +
    coord_cartesian(xlim = XLIM, expand = FALSE) +
    scale_y_discrete(limits = rev(levels(df_g$motif_num))) +
    theme_void() +
    theme(
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    ggtitle(paste0(g, n_text))
  
  out_png_heat <- file.path(
    ANALYSIS_DIR,
    paste0("motif_heatmap_binned_", safe_filename(g), "_", suffix, "_heat.png")
  )
  
  ggsave(
    out_png_heat,
    p_g_heat,
    width  = 6,
    height = 10,
    dpi    = 300
  )
  
  message("Saved minimal heatmap PNG for group: ", g, " -> ", basename(out_png_heat))
}

# Combined heatmap panel (all groups, full versions)
p_heat <- patchwork::wrap_plots(peth_plots_full, ncol = 4)
message("Heatmap generation done.")

# ------------------ PLOT: MEAN ± SEM PER MOTIF ------------------
SMOOTH_K <- 5  # window size in bins 

agg_motif_binned_smooth <- agg_motif_binned %>%
  dplyr::group_by(group, motif_num, signal) %>%
  dplyr::arrange(nor_bin, .by_group = TRUE) %>%
  dplyr::mutate(
    mean_smooth = zoo::rollapply(
      mean_value,
      width   = SMOOTH_K,
      FUN     = function(x) mean(x, na.rm = TRUE),
      fill    = NA,
      align   = "center"
    ),
    sem_smooth = zoo::rollapply(
      sem_value,
      width   = SMOOTH_K,
      FUN     = function(x) mean(x, na.rm = TRUE),
      fill    = NA,
      align   = "center"
    )
  ) %>%
  dplyr::ungroup()

group_colors <- c(
  "dLight_dStr" = "#8db4d6",
  "dLight_vStr" = "#83cdb6",
  "anxa1-flp"    = "#eb8c47",
  "dat-cre"      = "#249394"
)

p_motif_smooth <- ggplot(
  agg_motif_binned_smooth %>% dplyr::filter(signal == "signal"),
  aes(x = nor_bin, y = mean_smooth, colour = group, fill = group)
) +
  # SEM ribbon
  geom_ribbon(
    aes(
      ymin = mean_smooth - sem_smooth,
      ymax = mean_smooth + sem_smooth
    ),
    alpha  = 0.25,
    colour = NA
  ) +
  # Mean line
  geom_line(size = 0.9) +
  # Vertical line at 0
  geom_vline(xintercept = 0, linetype = "dashed") +
  # Facet: 6 columns
  facet_wrap(~ motif_num, ncol = 5) +
  # Colors
  scale_colour_manual(values = group_colors) +
  scale_fill_manual(values   = group_colors) +
  # Labels
  labs(
    x     = "Normalized time",
    y     = "Signal (mean ± SEM, smoothed)",
    title = "Motif-triggered signal (smoothed)"
  ) +
  theme_bw() +
  theme(legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90")
  )

p_motif_smooth_all <- ggplot(
  agg_motif_binned_smooth %>% dplyr::filter(signal == "signal" & group == "anxa1-flp"),
  aes(x = nor_bin, y = mean_smooth, colour = factor(motif_num))) +
  # Mean line
  geom_line(size = 0.9) +
  # Vertical line at 0
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # Labels
  labs(
    x     = "Normalized time",
    y     = "Signal (mean ± SEM, smoothed)",
    title = "Motif-triggered signal (smoothed)"
  ) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90"))

ggsave(
  file.path(ANALYSIS_DIR, "sup_figure_6F.pdf"),
  plot     = p_motif_smooth_all,
  width    = 6,
  height   = 4,
  units    = "in"
)

p_speed_smooth <- ggplot(
  agg_motif_binned_smooth %>% dplyr::filter(signal == "speed"),
  aes(x = nor_bin, y = mean_smooth, colour = group, fill = group)
) +
  # SEM ribbon
  geom_ribbon(
    aes(
      ymin = mean_smooth - sem_smooth,
      ymax = mean_smooth + sem_smooth
    ),
    alpha  = 0.25,
    colour = NA
  ) +
  # Mean line
  geom_line(size = 0.9) +
  # Vertical line at 0
  geom_vline(xintercept = 0, linetype = "dashed") +
  # Facet: 6 columns
  facet_wrap(~ motif_num, ncol = 5) +
  # Colors
  scale_colour_manual(values = group_colors) +
  scale_fill_manual(values   = group_colors) +
  # Labels
  labs(
    x     = "Normalized time",
    y     = "Speed (mean ± SEM, z-scored)",
    title = "Motif-triggered signal (smoothed)"
  ) +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90")
  )

sup_figure <- p_motif_smooth / p_speed_smooth

# Save as 12x12 inch PDF
ggsave(
  file.path(ANALYSIS_DIR, "sup_figure_lines.pdf"),
  plot     = sup_figure,
  width    = 12,
  height   = 24,
  units    = "in"
)

# ------------------ 5) AUC (0–1 s) PER GROUP × MOTIF ------------------
OUT_AUC <- file.path(ANALYSIS_DIR, "motif_auc_0to1_by_group.csv")

auc_long <- agg_motif_binned %>%
  dplyr::filter(signal == "signal") %>%
  dplyr::filter(nor_bin >= 0, nor_bin <= 1) %>%
  dplyr::group_by(group, motif_num, nor_bin) %>%
  dplyr::summarise(
    val              = mean(mean_value, na.rm = TRUE),
    mean_speed_motif = suppressWarnings(first(na.omit(mean_speed_motif))),
    .groups          = "drop_last"
  ) %>%
  dplyr::arrange(group, motif_num, nor_bin) %>%
  dplyr::summarise(
    mean_speed_motif = suppressWarnings(first(na.omit(mean_speed_motif))),
    auc = {
      v  <- val
      t  <- nor_bin
      dt <- dplyr::lead(t) - t
      trap <- (v + dplyr::lead(v)) / 2 * dt
      sum(trap, na.rm = TRUE)
    },
    .groups = "drop"
  ) %>%
  dplyr::rename(motif = motif_num)

auc_wide <- auc_long %>%
  tidyr::pivot_wider(
    id_cols    = c(motif, mean_speed_motif),
    names_from = group,
    values_from = auc,
    names_prefix = "auc_"
  ) %>%
  dplyr::arrange(dplyr::desc(mean_speed_motif), motif)

readr::write_csv(auc_wide, OUT_AUC)
message("Saved AUC table: ", OUT_AUC)

# ------------------ 6) SCATTER PLOTS & CORRELATIONS ------------------
# Standardize column name if needed
if ("auc_anxa1-flp" %in% names(auc_wide)) {
  auc_wide <- auc_wide %>%
    dplyr::rename(auc_anxa1_flp = `auc_anxa1-flp`)
}

# 1) Mean speed vs AUC (anxa1-flp)
p_anxa1_vs_meanspeed <- ggpubr::ggscatter(
  auc_wide,
  x   = "mean_speed_motif",
  y   = "auc_anxa1_flp",
  add = "reg.line", 
  conf.int = TRUE
) +
  ggpubr::stat_cor(
    method  = "pearson",
    label.x = min(auc_wide$mean_speed_motif, na.rm = TRUE),
    label.y = max(auc_wide$auc_anxa1_flp,      na.rm = TRUE)
  ) +
  labs(
    x = "Mean speed (motif)",
    y = "AUC — anxa1-flp"
  )

# 2) AUC(dLight dStr) vs AUC(anxa1-flp)
p_anxa1_vs_dStr <- ggpubr::ggscatter(
  auc_wide,
  x   = "auc_dLight_dStr",
  y   = "auc_anxa1_flp",
  add = "reg.line", 
  conf.int = TRUE
) +
  ggpubr::stat_cor(
    method  = "pearson",
    label.x = min(auc_wide$auc_dLight_dStr, na.rm = TRUE),
    label.y = max(auc_wide$auc_anxa1_flp,      na.rm = TRUE)
  ) +
  labs(
    x = "AUC — dLight_dStr",
    y = "AUC — anxa1-flp"
  )

# 3) AUC(dLight vStr) vs AUC(anxa1-flp)
p_anxa1_vs_vStr <- ggpubr::ggscatter(
  auc_wide,
  x   = "auc_dLight_vStr",
  y   = "auc_anxa1_flp",
  add = "reg.line", 
  conf.int = TRUE
) +
  ggpubr::stat_cor(
    method  = "pearson",
    label.x = min(auc_wide$auc_dLight_vStr, na.rm = TRUE),
    label.y = max(auc_wide$auc_anxa1_flp,      na.rm = TRUE)
  ) +
  labs(
    x = "AUC — dLight_vStr",
    y = "AUC — anxa1-flp"
  )

p_scatter <- p_anxa1_vs_meanspeed /
  p_anxa1_vs_dStr /
  p_anxa1_vs_vStr

# ------------------ 7) FINAL 4-PANEL PATCHWORK ------------------
p_all <- p_len | p_speed | p_heat | p_scatter
p_all <- p_all + patchwork::plot_layout(widths = c(0.5, 0.5, 4, 1))

ggsave(
  file.path(ANALYSIS_DIR, "figure_5DEF.pdf"),
  p_all,
  width  = 18,
  height = 10,
  dpi    = 300
)

message("All analyses done.")
