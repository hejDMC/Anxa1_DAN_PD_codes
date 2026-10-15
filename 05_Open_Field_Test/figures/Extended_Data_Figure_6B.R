# ============================================================
# Extended Data Figure 6B — Example aligned behavioral traces
# Author: Alessandro Contestabile
#
# Description:
#   - Visualizes an example open-field recording combining fiber-photometry activity and behavioral measurements.
#   - Displays synchronized photometry signal, locomotor speed, angular speed, and VAME motif identities.
#   - Represents consecutive behavioral motif episodes as a color-coded ethogram.
#   - Generates a four-panel figure illustrating the temporal alignment of neural and behavioral data.
#
# Requirements:
#   - Frame-wise *_aligned.csv files generated during OFT preprocessing.
#   - Columns: time, signal, speed, angular.speed, motif.
#
# Outputs:
#   - sup_figure_ethogram.pdf
#
# Notes:
#   - Uses a representative recording selected by video ID.
#   - Displays the interval between 352 and 368 seconds.
#   - Motif colors are assigned using the Spectral palette.
# ============================================================


# ------------------ LIBRARIES ------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(RColorBrewer)
})

# ------------------ SETTINGS ------------------
OUT_ROOT <- "path/to/VAME_project/results"

# Analysis directory for all outputs
ANALYSIS_DIR <- file.path(OUT_ROOT, "analysis")
dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------ PICK ONE ALIGNED FILE ------------------
aligned_files <- list.files(
  path      = OUT_ROOT,
  pattern   = "_aligned\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (!length(aligned_files)) {
  stop("No '*_aligned.csv' files found under OUT_ROOT: ", OUT_ROOT)
}

example_file <- aligned_files[grepl("tmaze_2023-05-16T12_21_10", aligned_files)][1]
message("Using example aligned file: ", example_file)

# ------------------ LOAD DATA ------------------
df <- readr::read_csv(example_file, show_col_types = FALSE)

needed_cols  <- c("time", "signal", "speed", "angular.speed", "motif")
missing_cols <- setdiff(needed_cols, names(df))
if (length(missing_cols)) {
  stop("Missing columns in example file: ", paste(missing_cols, collapse = ", "))
}

df <- df %>% arrange(time)

# ------------------ MOTIF EPISODES FOR COLOR STRIP ------------------
r   <- rle(df$motif)
len <- r$lengths
val <- r$values

end_idx   <- cumsum(len)
start_idx <- c(1L, head(end_idx, -1) + 1L)

motif_segments <- tibble(
  motif = factor(val),
  xmin  = df$time[start_idx],
  xmax  = df$time[end_idx]
)

# Spectral palette extended to number of motifs
n_motifs <- length(levels(motif_segments$motif))
base_pal <- brewer.pal(min(11, max(3, n_motifs)), "Spectral")
spectral_long <- colorRampPalette(base_pal)(n_motifs)
names(spectral_long) <- levels(motif_segments$motif)

# ------------------ PLOTS ------------------

# Signal (z-score)
p_signal <- ggplot(df, aes(x = time, y = signal)) +
  geom_line(linewidth = 0.4) +
  xlim(352, 368) +
  labs(
    y = "Signal (z)",
    x = NULL,
    title = "Fiber-photometry signal"
  ) +
  theme_classic()

# Speed (cm/s)
p_speed <- ggplot(df, aes(x = time, y = speed)) +
  xlim(352, 368) +
  geom_line(linewidth = 0.4) +
  labs(
    y = "Speed (cm/s)",
    x = NULL,
    title = "Locomotor speed"
  ) +
  theme_classic()

# Angular speed
p_ang <- ggplot(df, aes(x = time, y = angular.speed)) +
  xlim(352, 368) +
  geom_line(linewidth = 0.4) +
  labs(
    y = "Angular speed",
    x = NULL,
    title = "Angular speed (|Δθ|/s)"
  ) +
  theme_classic()

# Motif strip: colored blocks over time (no grids/axes, Spectral palette)
p_motif <- ggplot(motif_segments) +
  xlim(352, 368) +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = 0,
      ymax = 1,
      fill = motif
    ),
    color = NA
  ) +
  scale_fill_manual(values = spectral_long, guide = "none") +
  labs(x = "Time (s)", y = NULL, title = "Motifs") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, vjust = 0, size = 11),
    plot.margin = margin(t = 2, r = 5, b = 2, l = 5)
  )

# ------------------ COMBINE (4 ROWS) ------------------
p_example <- p_signal / p_speed / p_ang / p_motif +
  plot_layout(heights = c(1, 1, 1, 0.4)) +
  plot_annotation(title = "Example aligned trace: signal, kinematics, motifs")

# ------------------ SAVE EXAMPLE PLOT ------------------
out_pdf_example <- file.path(ANALYSIS_DIR, "sup_figure_ethogram.pdf")

ggsave(
  filename = out_pdf_example,
  plot     = p_example,
  width    = 12,      # adjust as you prefer
  height   = 6,       # adjust as you prefer
  dpi      = 300
)

message("Saved example plot: ", out_pdf_example)