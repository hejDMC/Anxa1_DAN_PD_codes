# MERSCOPE analysis

Python scripts for quality control, clustering, and visualization of
MERSCOPE spatial transcriptomics data from mouse ventral midbrain.

## Software requirements

The analysis uses Python, Scanpy, AnnData, NumPy, pandas, SciPy,
Matplotlib, and tifffile. Optional components require BBKNN and Zarr.

## Input data

The main pipeline requires an AnnData-compatible `.h5ad` or `.hdf5`
file containing:

- A cell-by-gene transcript-count matrix.
- Gene symbols in `var_names`.
- Unique cell identifiers in `obs_names`.
- Section identifiers in `obs["Datasets"]`.

Spatial overlay scripts additionally require cell coordinates, the
corresponding DAPI mosaic TIFF, and the appropriate coordinate transform.

Data availability: The MERSCOPE data have been deposited in the Gene Expression Omnibus (GEO) under accession number GSE351703: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE351703

## Scripts

- `01.Merscope_main_analysis.py`: quality control, normalization,
  whole-dataset clustering, TH/Slc6a3-based parent-cluster selection,
  and dopaminergic subclustering.
- `02.Plot_subcluster_marker_dotplots.py`: plots selected genes across subclusters.
- `03.Gene expression over dapi sections.py`: overlays gene expression on DAPI images.
- `04.Merscope_ overlay_clusters_on_dapi_sections.py`: overlays selected subclusters on DAPI images.
- `05.UMAP_per_condition_clusterofinterest.py`: plots existing subclusters
  separately by experimental group.
- `06.Figure_1H.py`: performs MASC-style mixed-effects analysis of dopaminergic neuron cluster proportions across MitoPark disease stages, reporting log2 odds ratios, 95% confidence intervals, and FDR-adjusted p-values for Figure 1H.

## Running the analysis

1. Update the input and output paths in each script.
2. Verify the section-to-animal mapping against the sample metadata.
3. Run `01.Merscope_main_analysis.py`. In Spyder, run the required
   cells in order; cells marked optional are exploratory or diagnostic.
4. Run the plotting scripts using the saved subcluster AnnData file.
5. For tissue overlays, select the matching section, DAPI image,
   and coordinate transform.

Cluster IDs must be checked after rerunning clustering; numerical labels
are not guaranteed to remain identical.

## Analysis overview

Blank probes are removed, and cells with fewer than 25 total counts
across retained genes are excluded. Counts are normalized per cell,
log-transformed, and scaled before PCA. A neighbor graph is used for
UMAP visualization and Leiden clustering.

A candidate dopaminergic parent cluster is selected using Th and
Slc6a3 expression and then subclustered.

The analysis uses the parameter values specified in the configuration
sections of the provided scripts.

## Expression plots

Gene-expression plots use the verified `log1p` layer, containing
unscaled, log-normalized expression.

For dotplots, dot size represents the fraction of cells with expression
above zero. When `standard_scale="var"` is enabled, colours are rescaled
separately for each gene and indicate relative expression patterns.

## Outputs

Scripts generate clustered AnnData files, cluster annotations, UMAPs,
marker dotplots, and spatial overlays. Output locations are configured
within each script.
