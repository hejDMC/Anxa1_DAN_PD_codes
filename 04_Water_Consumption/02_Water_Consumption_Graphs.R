# ============================================================
# Water Consumption — Fiber Photometry PETH and AUC Analysis
# ============================================================
#
# Description:
#   Analyze mouth-triggered fiber photometry responses during water consumption and quantify activity across time windows.
#
#   - Plot event-aligned PETHs (mean ± SEM)
#   - Calculate mean PETH curves for individual mice
#   - Compute per-mouse AUC using trapezoidal integration
#   - Compare AUCs across baseline, pre-event, and post-event
#     windows using paired statistical tests
#   - Combine PETH and AUC panels into a publication figure
#
# Input:
#   results/PETH_all.rds (or PETH_all.csv)
#
# Output:
#   results/figure_4/PETH_mouse_mean_curves.csv
#   results/figure_4/PETH_mouse_mean_AUCs.csv
#   results/figure_4/Figure_4DE.pdf
#
# Author: Alessandro Contestabile
# ============================================================

# ------------------ CLEAN ENVIRONMENT ------------------
rm(list = ls()); invisible(gc())

# ------------------ PATHS ------------------
wd <- "path/to/data"
stopifnot(dir.exists(wd))
setwd(wd)

res_dir  <- "results"
fig_dir  <- file.path(res_dir, "figure_4")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------ LIBRARIES ------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
})

# ------------------ PARAMETERS (FIGURE-SPECIFIC) ------------------
# AUC windows are defined on the PER-MOUSE mean curve (one curve per mouse × condition).
auc_windows <- tibble::tribble(
  ~window,     ~t0,  ~t1,
  "baseline",  -2,   -1,
  "-1_0",      -1,    0,
  "0_1",        0,    1
)

group_colors <- c(
  "dLight_dStr" = "#8db4d6",
  "dLight_vStr" = "#83cdb6",
  "fp_anxa1"    = "#eb8c47",
  "fp_dat"      = "#249394"
)

# ------------------ I/O HELPERS ------------------
load_peth <- function(base_name, dir = res_dir) {
  rds_path <- file.path(dir, paste0(base_name, ".rds"))
  csv_path <- file.path(dir, paste0(base_name, ".csv"))
  
  if (file.exists(rds_path)) return(readRDS(rds_path))
  if (file.exists(csv_path)) return(readr::read_csv(csv_path, show_col_types = FALSE))
  
  stop("Missing input: expected ", rds_path, " or ", csv_path)
}

make_safe <- function(df) {
  df <- df %>% mutate(across(where(is.factor), as.character))
  names(df) <- trimws(names(df))
  df
}

# ------------------ LOAD ------------------
PETH_all <- make_safe(load_peth("PETH_all"))

# ------------------ BASIC QC / TYPES ------------------
needed <- c("nor.time", "zscore", "event.type", "group", "mouse.ID", "body_part")
stopifnot(all(needed %in% names(PETH_all)))

PETH_all <- PETH_all %>%
  mutate(
    nor.time   = suppressWarnings(as.numeric(nor.time)),
    zscore     = suppressWarnings(as.numeric(zscore)),
    group      = as.character(group),
    event.type = as.character(event.type),
    mouse.ID   = as.character(mouse.ID),
    body_part  = as.character(body_part)
  )

# Consistent facet order (edit if you have other labels)
ev_levs <- c("rewarded", "empty", "error")
PETH_all <- PETH_all %>%
  mutate(event.type = factor(event.type, levels = intersect(ev_levs, unique(event.type)), ordered = TRUE))

# ------------------ SUMMARY FUNCTION (MEAN ± SEM) ------------------
mean_se <- function(x) {
  m  <- mean(x, na.rm = TRUE)
  sd <- stats::sd(x, na.rm = TRUE)
  n  <- sum(is.finite(x))
  se <- if (n > 0) sd / sqrt(n) else NA_real_
  c(y = m, ymin = m - se, ymax = m + se)
}

# ------------------ AUC HELPERS (TRAPEZOID) ------------------
clip_curve_to_window <- function(df_curve, t0, t1) {
  # df_curve: nor.time, mean_z
  df_curve <- df_curve %>% arrange(nor.time)
  
  seg <- df_curve %>% filter(nor.time >= t0, nor.time <= t1)
  
  # interpolate boundaries if not present
  if (!any(seg$nor.time == t0)) {
    y0  <- approx(df_curve$nor.time, df_curve$mean_z, xout = t0, ties = "ordered")$y
    seg <- bind_rows(tibble(nor.time = t0, mean_z = y0), seg)
  }
  if (!any(seg$nor.time == t1)) {
    y1  <- approx(df_curve$nor.time, df_curve$mean_z, xout = t1, ties = "ordered")$y
    seg <- bind_rows(seg, tibble(nor.time = t1, mean_z = y1))
  }
  
  seg %>% arrange(nor.time)
}

auc_from_segment <- function(seg) {
  if (nrow(seg) < 2) return(tibble(auc_raw = NA_real_))
  
  dt <- diff(seg$nor.time)
  y1 <- head(seg$mean_z, -1)
  y2 <- tail(seg$mean_z, -1)
  
  auc_raw <- sum(dt * (y1 + y2) / 2, na.rm = TRUE)
  tibble(auc_raw = auc_raw)
}

# ============================================================
# DATA FOR FIGURE: restrict to mouth + exclude event.type == error
# ============================================================
df_mouth <- PETH_all %>%
  filter(
    is.finite(nor.time), is.finite(zscore),
    body_part == "mouth",
    event.type != "error"
  )

# ------------------ PER-MOUSE MEAN CURVES ------------------
PETH_mouse_curve <- df_mouth %>%
  group_by(mouse.ID, group, event.type, nor.time) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

# ------------------ PER-MOUSE AUC ON MEAN CURVE ------------------
PETH_mouse_auc <- PETH_mouse_curve %>%
  group_by(mouse.ID, group, event.type) %>%
  group_modify(~{
    curve <- .x %>% arrange(nor.time)
    pmap_dfr(auc_windows, function(window, t0, t1) {
      seg  <- clip_curve_to_window(curve, t0, t1)
      auc  <- auc_from_segment(seg)
      bind_cols(tibble(window = window, t0 = t0, t1 = t1), auc)
    })
  }) %>%
  ungroup() %>%
  filter(is.finite(auc_raw)) %>%
  mutate(window = factor(window, levels = auc_windows$window))

# Save source data (publication-friendly)
write_csv(PETH_mouse_curve, file.path(fig_dir, "PETH_mouse_mean_curves.csv"))
write_csv(PETH_mouse_auc,   file.path(fig_dir, "PETH_mouse_mean_AUCs.csv"))

# ------------------ STATS SETUP (PAIRED WINDOW COMPARISONS) ------------------
# Pairing is within mouse.ID (same mouse contributes to all windows).
win_levels     <- levels(PETH_mouse_auc$window)
win_comparisons <- combn(win_levels, 2, simplify = FALSE)

# ============================================================
# PLOT — Panel A: PETH (episodes pooled)
# ============================================================
p_peth <- ggplot(
  df_mouth,
  aes(x = nor.time, y = zscore, colour = group, fill = group)
) +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.18, colour = NA) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  facet_grid(. ~ event.type) +
  scale_colour_manual(values = group_colors) +
  scale_fill_manual(values   = group_colors) +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(
    title = "Mouth-aligned PETH (episodes pooled)",
    x = "Time from event (s)",
    y = "z-score (ΔF/F0)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid     = element_blank(),
    strip.text     = element_text(face = "bold"),
    plot.title     = element_text(face = "bold"),
    legend.position = "bottom"
  )

# ============================================================
# PLOT — Panel B: AUC across windows (paired lines + paired tests)
# ============================================================
p_auc_win <- ggplot(
  PETH_mouse_auc,
  aes(x = window, y = auc_raw, colour = group, fill = group)
) +
  geom_boxplot(alpha = 0.20, outlier.shape = 16, outlier.size = 1.2) +
  
  # Paired lines: connect the same mouse across windows
  geom_line(
    aes(group = mouse.ID),
    alpha = 0.35,
    linewidth = 0.4
  ) +
  
  geom_point(
    size = 1.6, alpha = 0.75,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  
  facet_grid(group ~ event.type) +
  scale_colour_manual(values = group_colors) +
  scale_fill_manual(values   = group_colors) +
  
  # Paired stats within each facet (group × event.type)
  ggpubr::stat_compare_means(
    comparisons = win_comparisons,
    paired      = TRUE,
    method      = "t.test",       # switch to "wilcox.test" if preferred
    label       = "p.signif",
    hide.ns     = TRUE,
    size        = 3
  ) +
  
  labs(
    title = "Per-mouse AUC across windows (paired)",
    x = "Window",
    y = "AUC (mean curve)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    strip.text      = element_text(face = "bold"),
    plot.title      = element_text(face = "bold"),
    legend.position = "none"
  )

# ============================================================
# COMBINE (2:1 width) + SAVE
# ============================================================
combined <- p_peth + p_auc_win + plot_layout(widths = c(2, 1))

ggsave(
  filename = file.path(fig_dir, "Figure_4DE.pdf"),
  plot     = combined,
  width    = 12, height = 8, units = "in", dpi = 300
)

combined