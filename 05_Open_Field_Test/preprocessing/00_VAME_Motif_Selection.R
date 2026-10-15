
# ============================================================
# VAME motif frequency analysis
# Author: Alessandro Contestabile
#
# Description:
#   - Imports VAME motif labels from a 100-motif model.
#   - Calculates the frequency of each motif across recordings.
#   - Identifies motifs representing more than 1% of frames.
#   - Generates a motif frequency distribution plot.
#
# Requirements:
#   - VAME .npy motif labels generated using k = 100.
#
# Outputs:
#   - motif_frequency_100.csv
# ============================================================

# ------------------ CLEAN ENVIRONMENT ------------------
rm(list = ls())
invisible(gc())

# ------------------ LIBRARIES ------------------
suppressPackageStartupMessages({
  library(reticulate)
  library(ggplot2)
  library(dplyr)
})

# ------------------ SETTINGS ------------------
wd <- "path/to/VAME_project"

n.motif <- 100
frequency_threshold <- 1  # Percentage of total frames

output_dir <- file.path(wd, "motif_selection")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------ PYTHON / NUMPY ------------------
np <- reticulate::import("numpy", convert = TRUE)

# ------------------ IMPORT VAME MOTIFS ------------------
motif_files <- list.files(
  path = wd,
  pattern = paste0("^", n.motif, "_km_label_.*\\.npy$"),
  recursive = TRUE,
  full.names = TRUE
)

if (length(motif_files) == 0) {
  stop("No VAME motif files found.")
}

cat("Found", length(motif_files), "VAME files.\n")

motif_list <- lapply(seq_along(motif_files), function(i) {
  
  motifs <- as.integer(np$load(motif_files[i]))
  
  cat("[", i, "/", length(motif_files), "] ",
      basename(motif_files[i]), " imported\n", sep = "")
  
  motifs
})

all_motifs <- unlist(motif_list, use.names = FALSE)

# ------------------ MOTIF FREQUENCY ------------------
motif_frequency <- data.frame(
  motif = 0:(n.motif - 1),
  n_frames = tabulate(all_motifs + 1L, nbins = n.motif)
) %>%
  mutate(
    frequency_pct = 100 * n_frames / sum(n_frames),
    above_threshold = frequency_pct > frequency_threshold
  )

n_above_threshold <- sum(motif_frequency$above_threshold)

cat("\n========================================\n")
cat("VAME MOTIF FREQUENCY ANALYSIS\n")
cat("========================================\n")
cat("Total frames:", length(all_motifs), "\n")
cat("Number of motifs:", n.motif, "\n")
cat("Frequency threshold:", frequency_threshold, "%\n")
cat("Motifs above threshold:", n_above_threshold, "\n")

# ------------------ SAVE SUMMARY ------------------
write.csv(
  motif_frequency,
  file.path(output_dir, "motif_frequency_100.csv"),
  row.names = FALSE
)

# ------------------ PLOT MOTIF DISTRIBUTION ------------------
p <- ggplot(
  motif_frequency,
  aes(
    x = reorder(factor(motif), n_frames),
    y = n_frames
  )
) +
  geom_point(size = 2) +
  geom_hline(
    yintercept = length(all_motifs) * frequency_threshold / 100,
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    title = "VAME motif frequency distribution",
    subtitle = paste(
      n_above_threshold,
      "motifs above",
      frequency_threshold,
      "% of total frames"
    ),
    x = "Motif",
    y = "Number of frames"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

print(p)