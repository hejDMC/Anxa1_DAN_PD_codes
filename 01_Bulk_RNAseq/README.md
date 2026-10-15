# Bulk RNA-seq and single-nucleus RNA-seq analyses

This directory contains R scripts used to compare gene-expression changes in the MitoPark and conditional **Mfn2** knockout (cMfn2KO) models and to project the resulting bulk RNA-seq signatures onto midbrain dopaminergic nuclei from the publicly available **GSE233866** dataset (Salmani et al.). These analyses contribute to Figure 1 and Extended Data Figure 1 of the associated study.

## Contents

| File | Purpose |
| --- | --- |
| `01_MitoPark_MFN2_DEG_Overlap.R` | Compare precomputed bulk RNA-seq differential-expression results, identify genes downregulated in both models, and generate volcano, overlap, and fold-change plots. |
| `02_GSE233866_snRNAseq_Preprocessing.R` | Import and quality-control GSE233866 nuclei, perform Seurat normalization and clustering, select dopaminergic nuclei, and isolate the midbrain dopaminergic (mDA) region of interest. |
| `03_GSE233866_Bulk_Signature_Analysis.R` | Score mDA nuclei using the MitoPark, cMfn2KO, and shared downregulated-gene signatures; compare high- and low-score nuclei; generate plots and gene lists. |
| `mDA_selected_cells.txt` | Cell identifiers defining the manually selected mDA region of interest. Preserve this file to reproduce the original selection. |

## Input data

- **Bulk RNA-seq:** Precomputed differential-expression tables for the MitoPark and cMfn2KO models, supplied in `DEGs_mitopark_mfn2.xlsx` (worksheets `mitopark` and `mfn2`). Script 03 also expects `downregulated_DEGs_mitopark_mfn2.xlsx`, with the corresponding downregulated gene lists. The scripts analyze processed DEG tables; they do **not** perform alignment or differential-expression testing from raw sequencing reads.
- **GSE233866:** `GSE233866_metadata.csv` and `GSE233866_lesion_intact_counts.csv.gz`. Download the appropriate files from [NCBI GEO: GSE233866](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233866).
- **Intermediate object:** `mDA.rds`, produced by script 02 and used by script 03.

Input datasets and large Seurat `.rds` objects are not included in this directory.

## Running the analyses

Run the numbered scripts in order. Update the input/output directory variables near the beginning of each script to match your local installation. The preprocessing script uses `WD` as its working directory; the signature-analysis script also contains file paths that should be adjusted before execution.

1. **Bulk DEG overlap:** Run `01_MitoPark_MFN2_DEG_Overlap.R` with the bulk DEG workbook. It writes the combined plots and a CSV of shared downregulated genes.
2. **GSE233866 preprocessing:** Run `02_GSE233866_snRNAseq_Preprocessing.R` with the downloaded counts and metadata. It writes `allnuc.rds`, `DA.rds`, `mDA.rds`, and a nuclei QC summary.
3. **Signature projection:** Run `03_GSE233866_Bulk_Signature_Analysis.R` with `mDA.rds` and the two bulk DEG workbooks. It writes signature plots, individual gene-expression feature plots, and lists of available/missing genes.

### Reproducing the mDA selection

Script 02 originally defines the mDA region of interest using three interactive `CellSelector()` selections on a UMAP plot. The included `mDA_selected_cells.txt` records the resulting cell identifiers. For a non-interactive reproduction, replace the interactive selection block with:

```r
selected_cells <- readLines("path/to/mDA_selected_cells.txt")
selected_cells <- unique(selected_cells)
stopifnot(all(selected_cells %in% colnames(sobj)))
```


## Analysis notes

- The scripts retain the analysis-specific QC, clustering, module-score, and high/low score thresholds used for the study; see the scripts for exact parameters.
- The high/low groups are **defined by expression-derived module scores**. Differential-expression comparisons between these groups are descriptive and should not be interpreted as independent validation of the gene signatures.
- Interactive ROI selection and stochastic analysis steps can affect exact reproducibility. Preserve the cell-ID selection file and use consistent R/Seurat package versions and random seeds.

## Software

Analyses were performed in **R**, primarily using **Seurat**, **dplyr**, **ggplot2**, **readxl**, **data.table**, and associated plotting/normalization packages. Consult the `library()` statements in each script for the complete dependencies. Package versions are not pinned in this repository.
