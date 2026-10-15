
# ============================================================
# Extended Data Figure 8I — Locomotor speed distributions
# Author: Alessandro Contestabile
#
# Description:
#   - Imports frame-wise behavioral data from aligned open-field recordings.
#   - Compares locomotor speed distributions between control and experimental mice.
#   - Generates kernel density estimates of speed for each experimental cohort.
#   - Displays control distributions in grey and experimental distributions in cohort-specific colors.
#   - Arranges distributions by experiment for comparison.
#
# Requirements:
#   - Frame-wise *_aligned.csv files generated during OFT preprocessing.
#   - Experimental metadata and locomotor speed values.
#
# Outputs:
#   - Extended_Data_Figure_8I.pdf
#
# Notes:
#   - Only experiments specified in ALLOWED_EXPERIMENTS are included.
#   - Density estimates are calculated from pooled frame-wise speed measurements.
#   - The figure is descriptive; no statistical comparisons are performed in this script.
# ============================================================


# ------------------ CLEAN ENVIRONMENT ------------------
rm(list = ls())
invisible(gc())

# ------------------ LIBRARIES ------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

# ------------------ SETTINGS ------------------
wd <- "path/to/VAME_project/results"
setwd(wd)

OUTDIR <- file.path("analysis")
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

OUT_PDF <- file.path(OUTDIR, "Extended_Data_Figure_8I.pdf")

# Colors per experiment (used consistently across plots)
exp_colors <- c(
  "oft_6OHDA"               = "#F49495",
  "oft_mitopark_10-11weeks" = "#30C0C4",
  "oft_mitopark_15-18weeks" = "#0080C0",
  "oft_Tetx"                = "#ffb437",
  "oft_6ohda_tailStriatum" = "#ab3843ff",
  "oft_anti_anxa1"          = "#993896"  
)
ALLOWED_EXPERIMENTS <- names(exp_colors)

CTRL_COL <- "grey60"

# ------------------ HELPERS ------------------
safe_read_all_char <- function(f) {
  tryCatch(
    readr::read_csv(
      f,
      show_col_types = FALSE,
      progress = FALSE,
      col_types = readr::cols(.default = readr::col_character())
    ),
    error = function(e) {
      message("Failed reading: ", f, " | ", e$message)
      return(NULL)
    }
  )
}

pick_first_existing <- function(nms, candidates) {
  hit <- intersect(candidates, nms)
  if (!length(hit)) return(NA_character_)
  hit[1]
}

guess_ctrl_vs_exp <- function(df) {
  # If ctrl.vs.exp exists, normalize values
  if ("ctrl.vs.exp" %in% names(df)) {
    df$ctrl.vs.exp <- tolower(as.character(df$ctrl.vs.exp))
    df$ctrl.vs.exp <- ifelse(str_detect(df$ctrl.vs.exp, "ctrl|control|sham|wt|vehicle"), "ctrl", "exp")
    return(df)
  }
  
  # Otherwise infer from a common grouping/condition column
  candidate_cols <- intersect(
    c("group", "Group", "condition", "Condition", "genotype", "Genotype", "treatment", "Treatment"),
    names(df)
  )
  if (!length(candidate_cols)) {
    stop("No 'ctrl.vs.exp' column found and no grouping/condition/genotype column to infer it from.")
  }
  
  gcol <- candidate_cols[1]
  g <- tolower(as.character(df[[gcol]]))
  is_ctrl <- str_detect(g, "\\b(ctrl|control|sham|wt|vehicle)\\b")
  df$ctrl.vs.exp <- ifelse(is_ctrl, "ctrl", "exp")
  df
}

# ------------------ IMPORT ALIGNED FILES ------------------
aligned_files <- list.files(
  path = wd,
  pattern = "aligned\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(aligned_files)) stop("No files ending with 'aligned.csv' found under: ", wd)

dat_list <- lapply(aligned_files, safe_read_all_char)
dat_list <- dat_list[!vapply(dat_list, is.null, logical(1))]
if (!length(dat_list)) stop("All aligned files failed to read.")

dat <- dplyr::bind_rows(dat_list)

# ------------------ FIND REQUIRED COLUMNS ------------------
# experiment column
exp_col <- pick_first_existing(names(dat), c("experiment", "Experiment", "exp", "Exp"))
if (is.na(exp_col)) stop("No experiment column found. Expected one of: experiment/Experiment/exp/Exp.")

# speed/velocity column
speed_col <- pick_first_existing(names(dat), c("speed", "Speed", "velocity", "Velocity"))
if (is.na(speed_col)) stop("No speed/velocity column found. Expected one of: speed/Speed/velocity/Velocity.")

# Rename to standard names
dat <- dat %>%
  rename(
    experiment = all_of(exp_col),
    speed      = all_of(speed_col)
  )

# Ensure ctrl.vs.exp exists (or infer it)
dat <- guess_ctrl_vs_exp(dat)

# ------------------ FILTER + CLEAN ------------------
dat_filt <- dat %>%
  filter(experiment %in% ALLOWED_EXPERIMENTS) %>%
  mutate(
    experiment  = factor(experiment, levels = ALLOWED_EXPERIMENTS),
    ctrl.vs.exp = factor(ctrl.vs.exp, levels = c("ctrl", "exp")),
    speed_num   = suppressWarnings(as.numeric(speed))
  ) %>%
  filter(is.finite(speed_num))

if (!nrow(dat_filt)) stop("After filtering, no rows remain. Check experiment names and speed values.")

# ------------------ PLOT: DENSITY DISTRIBUTIONS ------------------
dat_filt <- dat_filt %>%
  mutate(color_key = ifelse(ctrl.vs.exp == "ctrl", "ctrl", as.character(experiment)))

col_map <- c("ctrl" = CTRL_COL, exp_colors)

p <- ggplot(dat_filt, aes(x = speed_num, color = color_key)) +
  geom_density(linewidth = 1.0, adjust = 1, alpha = 0.5) +
  facet_wrap(~ experiment, scales = "free_y", nrow = 1) +
  scale_color_manual(values = col_map, breaks = c("ctrl", ALLOWED_EXPERIMENTS)) +
  labs(
    x = "Velocity (speed)",
    y = "Density",
    title = "Full velocity distributions (CTRL vs EXP)",
    color = NULL
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

print(p)

# ------------------ SAVE ------------------
ggsave(filename = OUT_PDF, plot = p, width = 10, height = 4, units = "in", dpi = 300)
message("Saved PDF: ", OUT_PDF)
