# ============================================================
# GSE233866 — snRNA-seq Preprocessing and mDA Selection
# ============================================================
# Description:Preprocessing of the single-nucleus RNA-seq dataset GSE233866 (Salmani et al.).
#
#   - Import raw counts and sample metadata
#   - Perform nuclei quality control
#   - Normalize and cluster all nuclei
#   - Select dopaminergic nuclei (Slc6a3 or Th expression)
#   - Reprocess the dopaminergic population
#   - Select midbrain dopaminergic (mDA) nuclei by UMAP ROI
#   - Reprocess and save the final mDA population
#
# Input:
#   GSE233866_metadata.csv
#   GSE233866_lesion_intact_counts.csv.gz
#
# Output:
#   summary_stats_all_nuclei.csv
#   allnuc.rds
#   DA.rds
#   mDA_selected_cells.txt
#   mDA.rds
#
# Dataset: GEO GSE233866
# Author: Alessandro Contestabile
# ============================================================

# ------------------ CLEAN ENVIRONMENT ------------------
rm(list = ls()); invisible(gc())

# ------------------ LIBRARIES ------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(stringr)
  library(sctransform)
  library(future)
  library(scales)
  library(RColorBrewer)
  library(readxl)
  library(dplyr)
  library(dendextend)
  library(data.table)
  library(glmGamPoi)
  library(matrixStats)
})

# ------------------ HELPERS ------------------
rm_assay <- function(obj, assay = "SCT") {
  if (!(assay %in% Assays(obj))) return(obj)
  has_remove <- "SeuratObject" %in% loadedNamespaces() &&
    ("RemoveAssay" %in% getNamespaceExports("SeuratObject"))
  if (has_remove) {
    obj <- SeuratObject::RemoveAssay(obj, assay = assay)
  } else {
    obj[[assay]] <- NULL
  }
  obj
}

rm_clusters_graphs_reductions <- function(obj) {
  drop_cols <- grep("(^.*_snn_res\\.|^snn$|^seurat_clusters$)",
                    colnames(obj@meta.data), value = TRUE)
  if (length(drop_cols)) obj@meta.data[, drop_cols] <- NULL
  if (length(obj@graphs)) obj@graphs <- list()
  if (length(obj@reductions)) obj@reductions <- list()
  obj
}

# Build neighbors + clusters with a consistent graph name
make_clusters <- function(obj, dims = 1:30, k = 20, resolution = 0.4) {
  if ("SCT" %in% Assays(obj)) DefaultAssay(obj) <- "SCT"
  if (!"pca" %in% names(obj@reductions)) obj <- RunPCA(obj, npcs = max(dims), verbose = FALSE)
  if (!"umap" %in% names(obj@reductions)) obj <- RunUMAP(obj, dims = dims)
  gname <- if ("SCT" %in% Assays(obj)) "SCT_snn" else paste0(DefaultAssay(obj), "_snn")
  obj <- FindNeighbors(obj, reduction = "pca", dims = dims, k.param = k, graph.name = gname)
  obj <- FindClusters(obj, graph.name = gname, resolution = resolution, n.start = 100, n.iter = 100)
  Idents(obj) <- "seurat_clusters"
  obj
}

# ------------------ PARALLEL & REPRODUCIBILITY ------------------
if (.Platform$OS.type == "windows") {
  plan("multisession", workers = 4)
} else {
  plan("multicore", workers = 4)
}
options(future.globals.maxSize = 60 * 1024^3)  # 60 GB
set.seed(2301)

# ------------------ PATHS ------------------
WD <- "path/to/GSE233866"
setwd(WD)

# ------------------ LOAD METADATA ------------------
metadata <- read.csv("GSE233866_metadata.csv", stringsAsFactors = FALSE)
required_meta_cols <- c("Sample_Name","condition","age","gender","X6OHDA.dose.ug.ul")
stopifnot(all(required_meta_cols %in% colnames(metadata)))

meta <- metadata[, required_meta_cols]
meta$Sample_Name <- as.character(meta$Sample_Name)
meta$condition   <- as.character(meta$condition)
meta$age         <- as.character(meta$age)
meta$gender      <- as.character(meta$gender)

# ------------------ LOAD COUNTS and BASIC FILTER ------------------
GSE233866_1 <- fread("GSE233866_lesion_intact_counts.csv.gz", data.table = FALSE)
stopifnot(ncol(GSE233866_1) > 1L)

rownames(GSE233866_1) <- GSE233866_1[[1]]
GSE233866_1 <- GSE233866_1[, -1, drop = FALSE]
GSE233866_1 <- GSE233866_1[rowSums(GSE233866_1 > 0) >= 10, , drop = FALSE]

# ------------------ CREATE SEURAT OBJECT ------------------
sobj <- CreateSeuratObject(counts = GSE233866_1, project = "GSE233866_1")

# ------------------ ALIGN METADATA TO CELLS ------------------
sample_ids <- substr(colnames(GSE233866_1), 1, 4)
if (anyDuplicated(meta$Sample_Name)) {
  warning("Duplicated Sample_Name entries in metadata; proceeding by first match per sample ID.")
}
meta_aligned <- meta[match(sample_ids, meta$Sample_Name), , drop = FALSE]
if (anyNA(meta_aligned$Sample_Name)) {
  bad <- unique(sample_ids[is.na(meta_aligned$Sample_Name)])
  stop("Unmatched sample IDs (first 4 chars of cell names) not found in metadata$Sample_Name: ",
       paste(bad, collapse = ", "))
}
rownames(meta_aligned) <- colnames(sobj)
sobj <- AddMetaData(sobj, metadata = meta_aligned)
sobj$sample <- "GSE233866_1"

# ------------------ QC METRICS and BASIC CELL FILTERS ------------------
sobj[["percent.mt"]]   <- PercentageFeatureSet(sobj, pattern = "^mt-", assay = "RNA")
sobj[["percent.ribo"]] <- PercentageFeatureSet(sobj, pattern = "^(Rpl|Rps)", assay = "RNA")

sobj <- subset(
  sobj, subset = nFeature_RNA >= 500 & nFeature_RNA <= 10000 & percent.mt <= 5
)

if ("Malat1" %in% Features(sobj)) {
  sobj <- subset(sobj, features = setdiff(Features(sobj), "Malat1"))
}

message("Cells after QC: ", ncol(sobj), " (in Salmani et al.: 68,914)")

# ------------------ NORMALIZATION (SCTransform) ------------------
sobj <- SCTransform(sobj, variable.features.n = 1000)
Assays(sobj)

# ------------------ DIMENSION REDUCTION ------------------
DefaultAssay(sobj) <- "SCT"
sobj <- RunPCA(sobj, npcs = 50, verbose = FALSE)
sobj <- RunUMAP(sobj, dims = 1:30)

# ------------------ CLUSTERING (ALL NUCLEI) ------------------
sobj <- FindNeighbors(sobj, reduction = "pca", k.param = 20, dims = 1:30)
sobj <- FindClusters(sobj, resolution = seq(0.1, 0.9, by = 0.1), n.start = 100, n.iter = 100)
sobj <- SetIdent(sobj, value = "SCT_snn_res.0.1")

# ------------------ RNA LOG-NORM (FOR MARKERS/EXPORTS) ------------------
DefaultAssay(sobj) <- "RNA"
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)

# ------------------ QC SUMMARY (ALL NUCLEI) ------------------
umi_stats  <- colSums(GetAssayData(sobj, assay = "RNA", layer = "counts"))
gene_stats <- sobj$nFeature_RNA
expr_stats <- colSums(GetAssayData(sobj, assay = "RNA", layer = "data"))

umi_per_gene <- umi_stats / pmax(gene_stats, 1L)

summary_stats <- data.frame(
  total_cells    = ncol(sobj),
  umi_min        = as.integer(min(umi_stats)),
  umi_max        = as.integer(max(umi_stats)),
  umi_mean       = as.integer(round(mean(umi_stats))),
  genes_min      = as.integer(min(gene_stats)),
  genes_max      = as.integer(max(gene_stats)),
  genes_mean     = as.integer(round(mean(gene_stats))),
  umi_gene_mean  = round(mean(umi_per_gene), 2),
  umi_gene_range = paste0(round(range(umi_per_gene), 2), collapse = "–"),
  stringsAsFactors = FALSE
)

print(summary_stats)
write.csv(summary_stats, file = file.path(WD, "summary_stats_all_nuclei.csv"), row.names = FALSE)

# ------------------ SAVE ------------------
saveRDS(sobj, file = file.path(WD, "allnuc.rds"))

# ------------------ SETUP ------------------
DefaultAssay(sobj) <- "RNA"

# ------------------ GATE DA NUCLEI (SLC6A3|TH > 0, CASE-ROBUST) ------------------
genes_present <- rownames(sobj)
get_gene <- function(sym) {
  if (sym %in% genes_present) return(sym)
  alt <- genes_present[tolower(genes_present) == tolower(sym)]
  if (length(alt)) return(alt[1])
  NA_character_
}
g_Slc6a3 <- get_gene("Slc6a3")
g_Th     <- get_gene("Th")
if (anyNA(c(g_Slc6a3, g_Th))) stop("Required markers not found (Slc6a3/Th).")

rna_counts <- GetAssayData(sobj, assay = "RNA", layer = "counts")
da_mat <- rna_counts[c(g_Slc6a3, g_Th), , drop = FALSE]
da_cells <- colnames(da_mat)[colSums(da_mat > 0) > 0]

sobj <- subset(sobj, cells = da_cells)
message("DA-gated cells: ", ncol(sobj), " (in Salmani et al.: 33,052)")

# ------------------ CLEAR OLD NORMALIZATION/CLUSTERS ------------------
rm_assay <- function(obj, assay = "SCT") {
  if (!(assay %in% Assays(obj))) return(obj)
  has_remove <- "SeuratObject" %in% loadedNamespaces() &&
    ("RemoveAssay" %in% getNamespaceExports("SeuratObject"))
  if (has_remove) {
    obj <- SeuratObject::RemoveAssay(obj, assay = assay)
  } else {
    obj[[assay]] <- NULL
  }
  obj
}

sobj <- rm_assay(sobj, "SCT")
drop_cols <- grep("(^.*_snn_res\\.|^snn$|^seurat_clusters$)", colnames(sobj@meta.data), value = TRUE)
if (length(drop_cols)) sobj@meta.data[, drop_cols] <- NULL
if (length(sobj@graphs))     sobj@graphs     <- list()
if (length(sobj@reductions)) sobj@reductions <- list()

# ------------------ NORMALIZATION & DR (DA-ONLY) ------------------
set.seed(2301)
sobj <- SCTransform(sobj, variable.features.n = 1000)
DefaultAssay(sobj) <- "SCT"
sobj <- RunPCA(sobj, npcs = 50, verbose = FALSE)
sobj <- RunUMAP(sobj, dims = 1:30)

# ------------------ CLUSTERING (DA-ONLY) ------------------
sobj <- make_clusters(sobj, dims = 1:30, k = 20, resolution = 0.4)

# ------------------ RNA LOG-NORM  ------------------
DefaultAssay(sobj) <- "RNA"
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# ------------------ SAVE ------------------
saveRDS(sobj, file = file.path(WD, "DA.rds"))

# ------------------ SELECT REGION OF INTEREST ------------------
FeaturePlot(sobj, features = c("Ddc", "En1", "Nr4a2", "Prlr", "Slc17a7", "Slc17a6"), reduction = "umap", ncol = 3)
p <- FeaturePlot(sobj, features = g_Th, reduction = "umap")
selected_1 <- CellSelector(plot = p)
selected_2 <- CellSelector(plot = p)
selected_3 <- CellSelector(plot = p)
selected_cells <- unique(c(selected_1, selected_2, selected_3))
message("Selected cells in ROI: ", length(selected_cells))

sel_path <- file.path(WD, "mDA_selected_cells.txt")
writeLines(selected_cells, sel_path)

# ------------------ SUBSET TO ROI (mDA) ------------------
sobj <- subset(sobj, cells = selected_cells)
message("DA-gated cells: ", ncol(sobj), " (in Salmani et al.: 24,809)")

# ------------------ CLEAR OLD NORMALIZATION/CLUSTERS (mDA) ------------------
sobj <- rm_assay(sobj, "SCT")
drop_cols2 <- grep("(^.*_snn_res\\.|^snn$|^seurat_clusters$)", colnames(sobj@meta.data), value = TRUE)
if (length(drop_cols2)) sobj@meta.data[, drop_cols2] <- NULL
if (length(sobj@graphs))     sobj@graphs     <- list()
if (length(sobj@reductions)) sobj@reductions <- list()

# ------------------ NORMALIZATION & DR (mDA ROI) ------------------
sobj <- SCTransform(sobj, variable.features.n = 1000)
DefaultAssay(sobj) <- "SCT"
sobj <- RunPCA(sobj, npcs = 50, verbose = FALSE)
sobj <- RunUMAP(sobj, dims = 1:30)

# ------------------ CLUSTERING (mDA ROI) ------------------
sobj <- make_clusters(sobj, dims = 1:30, k = 20, resolution = 0.4)

# ------------------ RNA LOG-NORM (mDA FOR MARKERS/EXPORTS) ------------------
DefaultAssay(sobj) <- "RNA"
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# ------------------ SAVE ------------------
saveRDS(sobj, file = file.path(WD, "mDA.rds"))
