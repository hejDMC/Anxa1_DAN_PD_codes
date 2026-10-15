# ============================================================
# VAME motif similarity analysis and hierarchical clustering
# Author: Alessandro Contestabile
#
# Description:
#   - Imports VAME motif labels and latent representations
#     from individual Open Field Test (OFT) recordings.
#   - Computes the mean latent representation of each motif
#     across all recordings.
#   - Performs hierarchical clustering to assess similarities
#     between behavioral motifs.
#   - Applies predefined motif merges based on inspection
#     of motif similarity.
#   - Recomputes the dendrogram after motif merging.
#
# Requirements:
#   - VAME model outputs generated with k = 42 motifs.
#   - details.csv containing video identifiers.
#   - NumPy .npy files containing motif labels and
#     latent representations for each recording.
#
# Outputs:
#   - dendrogram_vame.pdf:
#     Hierarchical clustering of motifs after merging.
#
# Notes:
#   - Motif merging is manually specified and is not
#     performed automatically by the clustering algorithm.
#   - Hierarchical clustering uses Euclidean distances
#     between mean motif latent representations and
#     complete-linkage clustering (R hclust defaults).
#   - The original VAME motif assignments are not modified
#     on disk.
# ============================================================

# ------------------ CLEAN ENVIRONMENT ------------------
rm(list = ls())
invisible(gc())

# ------------------ LIBRARIES ------------------
suppressPackageStartupMessages({
  library(reticulate)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(rhdf5)
})

# ------------------ SETTINGS ------------------
frame.behavior <- 1/30
n.motif        <- 42

# Base working directory for VAME model
wd <- "path/to/VAME_project"
setwd(wd)

# ------------------ IMPORT DETAILS ------------------
details <- read.csv("details.csv")

# Select videos listed in the experimental metadata
videos  <- details$video

# ------------------ PYTHON / NUMPY HANDLE ------------------
np <- import("numpy")

# ------------------ FILE LIST & FILTER ------------------
# All .mp4 in videos/ folder
files.to.read <- list.files(
  path      = "videos",
  pattern   = "\\.mp4$",
  recursive = TRUE,
  full.names = TRUE
)

# Extract video names without extension
video_names <- sub("\\.mp4$", "", basename(files.to.read))

# Keep only those that are present in `details$video`
keep_idx        <- video_names %in% videos
files.to.read   <- files.to.read[keep_idx]
video_names     <- video_names[keep_idx]

cat("Found", length(files.to.read), "video(s) with matching metadata.\n")

# ------------------ IMPORT VAME DATA (FAST, NO RBIND IN LOOP) ------------------

# We store each video result in a list and bind at the end.
vame_list <- vector("list", length(files.to.read))

for (i in seq_along(files.to.read)) {
  video_path <- files.to.read[i]
  video      <- video_names[i]
  
  # Paths to .npy files
  base_res_dir <- file.path(
    wd, "results", video,
    paste0("VAME_", n.motif),
    paste0("hmm-",  n.motif)
  )
  
  motifs_file  <- file.path(base_res_dir, paste0(n.motif, "_km_label_",  video, ".npy"))
  latent_file  <- file.path(base_res_dir, paste0("latent_vector_",      video, ".npy"))
  
  if (!file.exists(motifs_file)) {
    warning("Missing motifs file for video: ", video, " (", motifs_file, ")")
    next
  }
  if (!file.exists(latent_file)) {
    warning("Missing latent_vector file for video: ", video, " (", latent_file, ")")
    next
  }
  
  # Load numpy arrays
  motifs_np  <- np$load(motifs_file)
  latent_np  <- np$load(latent_file)
  
  # Convert to data.frame
  vame_motifs <- data.frame(motifs = as.integer(motifs_np))
  latent_df   <- as.data.frame(latent_np)
  
  # Build final data.frame for this video
  vame_df <- cbind(
    time  = seq_len(nrow(vame_motifs)) - 1L,
    motif = vame_motifs$motifs,
    latent_df
  )
  
  vame_df$time  <- vame_df$time * frame.behavior
  vame_df$video <- video
  
  vame_list[[i]] <- vame_df
  
  cat("[", i, "/", length(files.to.read), "] ", video, " imported\n", sep = "")
}

# Remove NULL entries (in case of missing files) and bind
vame_list <- vame_list[!vapply(vame_list, is.null, logical(1))]
if (length(vame_list) == 0L) {
  stop("No VAME data imported. Check paths / filenames / details.")
}

vame_features <- dplyr::bind_rows(vame_list)
cat("Total rows in vame_features:", nrow(vame_features), "\n")

# ------------------ DENDROGRAM: BEFORE MERGING ------------------
latent_cols <- 3:32

vame_features.agg <- aggregate(
  vame_features[, latent_cols],
  by = list(motif = vame_features$motif),
  FUN = mean
)

rownames(vame_features.agg) <- vame_features.agg$motif

motifs_hc <- hclust(dist(vame_features.agg[, -1]))
plot(motifs_hc, hang = -1, cex = 1)

# ------------------ MOTIF MERGING / REMAPPING ------------------
# Predefined merges of motifs with similar latent representations.
# These mappings were selected manually following inspection of motif similarity and are applied consistently in downstream analyses.
vame_features$motif[vame_features$motif == 11] <- 28
vame_features$motif[vame_features$motif == 2] <- 8
vame_features$motif[vame_features$motif == 15] <- 38
vame_features$motif[vame_features$motif == 29] <- 34

# ------------------ DENDROGRAM: AFTER MERGING ------------------
vame_features.agg2 <- aggregate(
  vame_features[, latent_cols],
  by = list(motif = vame_features$motif),
  FUN = mean
)

rownames(vame_features.agg2) <- vame_features.agg2$motif
motifs_hc2 <- hclust(dist(vame_features.agg2[, -1]))

pdf(file = file.path(wd, "dendrogram_vame.pdf"), width = 12, height = 6)
plot(motifs_hc2, hang = -1, cex = 1)
dev.off()

plot(motifs_hc2, hang = -1, cex = 1)
cat("Dendrogram saved to:", file.path(wd, "dendrogram_vame.pdf"), "\n")