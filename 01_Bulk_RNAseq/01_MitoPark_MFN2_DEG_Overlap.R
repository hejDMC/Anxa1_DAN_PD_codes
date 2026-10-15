# ============================================================
# MitoPark and cMfn2KO — Differential expression comparison
# Author: Alessandro Contestabile
#
# Description:
#   - Imports differential expression results from MitoPark and cMfn2KO experimental models.
#   - Identifies significantly downregulated genes in each dataset (adjusted p-value < 0.05).
#   - Generates volcano plots highlighting the most significant genes and Anxa1.
#   - Identifies genes downregulated in both models.
#   - Compares log2 fold changes of overlapping genes.
#   - Generates a combined figure containing volcano plots, a Venn diagram, and a fold-change scatter plot.
#
# Requirements:
#   - DEGs_mitopark_mfn2.xlsx
#   - Worksheets: "mitopark" and "mfn2"
#   - Columns including Symbol, log2FoldChange, and padj.
#
# Outputs:
#   - graphs_bulk_down.pdf
#
# Notes:
#   - Mitochondrial and selected ribosomal genes are excluded.
#   - Only genes with log2FoldChange <= 0 are retained for visualization.
#   - Significantly downregulated genes are defined by padj < 0.05 and log2FoldChange < 0.
#   - Overlapping genes are ranked by the sum of absolute log2 fold changes across both models.
# ============================================================

# ============================================================
# 1. SETUP
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggvenn)
  library(cowplot)
  library(ggrepel)
})

# ============================================================
# 2. PATHS AND PARAMETERS
# ============================================================

DATA_DIR <- "path/to/Bulk_RNAseq/data"
OUT_DIR  <- "path/to/Bulk_RNAseq/results"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DEG_FILE <- file.path(DATA_DIR, "DEGs_mitopark_mfn2.xlsx")
OUT_PDF  <- file.path(OUT_DIR, "graphs_bulk_down.pdf")
OUT_CSV  <- file.path(OUT_DIR, "MitoPark_MFN2_shared_downregulated_genes.csv")

PADJ_THRESHOLD <- 0.05
N_LABELS <- 10

MITOPARK_COLOR <- "#0080C0"
MFN2_COLOR     <- "#993292"
NS_COLOR       <- "lightgrey"


# ============================================================
# 3. LOAD AND CLEAN MITOPARK DATASET
# ============================================================
mitopark <- read_excel(DEG_FILE, sheet = "mitopark")
mitopark <- mitopark[, 1:8]
mitopark <- na.omit(mitopark)
colnames(mitopark)[1] <- "Symbol"

mitopark <- mitopark %>%
  filter(!grepl("^mt-|^rp[0-9]", Symbol, ignore.case = TRUE), log2FoldChange <= 0) %>%
  mutate(sig = ifelse(padj < PADJ_THRESHOLD, "sig", "ns"))

# Label the 10 most significant genes and Anxa1 (if significant)
mitopark_labeled <- mitopark %>%
  filter(sig == "sig") %>%
  arrange(padj) %>%
  slice_head(n = N_LABELS) %>%
  bind_rows(mitopark %>% filter(Symbol == "Anxa1", sig == "sig")) %>%
  distinct(Symbol, .keep_all = TRUE)

# ============================================================
# 4. MITOPARK VOLCANO PLOT
# ============================================================
p_volcano_mito <- ggplot(mitopark, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), alpha = 0.8, size = 1.5) +
  geom_text_repel(data = mitopark_labeled, aes(label = Symbol), size = 3, max.overlaps = 100) +
  scale_color_manual(values = c("sig" = MITOPARK_COLOR, "ns" = NS_COLOR)) +
  ylim(0, 45) +
  xlim(-4.5, 0) +
  geom_hline(yintercept = -log10(PADJ_THRESHOLD), linetype = "dashed", color = "darkred") +
  labs(title = "Mitopark Volcano Plot", x = "log2 Fold Change", y = "-log10(padj)") +
  theme_minimal() +
  theme(legend.position = "none")

# ============================================================
# 5. LOAD AND CLEAN MFN2 DATASET
# ============================================================

mfn2 <- read_excel(DEG_FILE, sheet = "mfn2")
mfn2 <- mfn2[, 1:8]
mfn2 <- na.omit(mfn2)
colnames(mfn2)[1] <- "Symbol"

mfn2 <- mfn2 %>%
  filter(!grepl("^mt-|^rp[0-9]", Symbol, ignore.case = TRUE), log2FoldChange <= 0) %>%
  mutate(sig = ifelse(padj < PADJ_THRESHOLD, "sig", "ns"))

# Label the 10 most significant genes and Anxa1 (if significant)
mfn2_labeled <- mfn2 %>%
  filter(sig == "sig") %>%
  arrange(padj) %>%
  slice_head(n = N_LABELS) %>%
  bind_rows(mfn2 %>% filter(Symbol == "Anxa1", sig == "sig")) %>%
  distinct(Symbol, .keep_all = TRUE)

# ============================================================
# 6. MFN2 VOLCANO PLOT
# ============================================================
p_volcano_mfn2 <- ggplot(mfn2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), alpha = 0.8, size = 1.5) +
  geom_text_repel(data = mfn2_labeled, aes(label = Symbol), size = 3, max.overlaps = 100) +
  scale_color_manual(values = c("sig" = MFN2_COLOR, "ns" = NS_COLOR)) +
  ylim(0, 45) +
  xlim(-4.5, 0) +
  geom_hline(yintercept = -log10(PADJ_THRESHOLD), linetype = "dashed", color = "darkred") +
  labs(title = "MFN2 Volcano Plot", x = "log2 Fold Change", y = "-log10(padj)") +
  theme_minimal() +
  theme(legend.position = "none")

# ============================================================
# 7. IDENTIFY SIGNIFICANTLY DOWNREGULATED GENES
# ============================================================
mitopark_sig <- mitopark %>%
  filter(padj < PADJ_THRESHOLD, log2FoldChange < 0)

mfn2_sig <- mfn2 %>%
  filter(padj < PADJ_THRESHOLD, log2FoldChange < 0)

cat("\nSignificantly downregulated genes:\n")
cat("MitoPark:", n_distinct(mitopark_sig$Symbol), "\n")
cat("MFN2:", n_distinct(mfn2_sig$Symbol), "\n")

# ============================================================
# 8. VENN DIAGRAM
# ============================================================
gene_sets <- list(
  Mitopark = mitopark_sig$Symbol,
  MFN2 = mfn2_sig$Symbol
)

p1 <- ggvenn(gene_sets,
             fill_color = c(MITOPARK_COLOR, MFN2_COLOR),
             stroke_size = 0.5,
             set_name_size = 5,
             show_percentage = FALSE) +
  coord_flip() +
  ggtitle("Overlap of Significantly Downregulated Genes")

# ============================================================
# 9. IDENTIFY AND RANK OVERLAPPING GENES
# ============================================================
overlap_genes <- intersect(mitopark_sig$Symbol, mfn2_sig$Symbol)

overlap_df <- mitopark_sig %>%
  filter(Symbol %in% overlap_genes) %>%
  rename(logFC_mitopark = log2FoldChange) %>%
  inner_join(
    mfn2_sig %>%
      filter(Symbol %in% overlap_genes) %>%
      rename(logFC_mfn2 = log2FoldChange),
    by = "Symbol"
  ) %>%
  mutate(combined_magnitude = abs(logFC_mitopark) + abs(logFC_mfn2)) %>%
  arrange(desc(combined_magnitude))

cat("\nShared significantly downregulated genes:", length(overlap_genes), "\n")
print(head(overlap_df, 10))

# Export shared genes for downstream analyses
write.csv(overlap_df, OUT_CSV, row.names = FALSE)

# ============================================================
# 10. SCATTER PLOT — LOG2 FOLD CHANGE COMPARISON
# ============================================================

p2 <- ggplot(overlap_df, aes(x = logFC_mitopark, y = logFC_mfn2, label = Symbol)) +
  geom_point(color = "#4D5BA6", size = 2.5) +
  geom_text(size = 4, check_overlap = FALSE, vjust = 1.5) +
  geom_abline(linetype = "dashed", color = "gray") +
  ylim(-1.75, -0.2) +
  xlim(-1.75, -0.2) +
  theme_minimal() +
  labs(
    title = paste0("Downregulated Genes Overlap (n = ", nrow(overlap_df), ")"),
    x = "log2FC (Mitopark)",
    y = "log2FC (MFN2)"
  )

# ============================================================
# 11. COMBINE AND SAVE FIGURE
# ============================================================

combo_plot <- plot_grid(
  p_volcano_mito, p_volcano_mfn2, p1, p2,
  ncol = 4,
  rel_widths = c(2, 2, 1, 2)
)

print(combo_plot)

ggsave(OUT_PDF, combo_plot, width = 18, height = 6)

cat("\nAnalysis completed.\n")
cat("Figure saved to:", OUT_PDF, "\n")
cat("Shared genes saved to:", OUT_CSV, "\n")