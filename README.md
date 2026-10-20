# Anxa1+ Dopamine Neurons — Analysis Code

This repository contains data, analysis code, and figures related to the experiments on Anxa1+ dopamine neurons (DANs) of the following publication:

Mantas, I., Contestabile, A., Skara, V., Loiseau, C., Rubino, E., Santos, I. A., Cramb, K. M. L., Filograna, R., Wade-Martins, R., Magill, P. J., & Meletis, K.  
**Anxa1+ dopamine neuron vulnerability defines prodromal Parkinson's disease-like bradykinesia and procedural motor learning impairment**

The analysis code is organized by experimental approach.

## Repository structure

- **`01_Bulk_RNAseq/`** — Bulk RNA-seq processing and analysis of gene-expression changes.
- **`02_MERFISH/`** — MERSCOPE spatial transcriptomics analysis, including quality control, clustering, marker visualization, spatial overlays, and cluster-abundance comparisons.
- **`03_Fiber_Photometry_Optogenetics/`** — Analysis and visualization of fiber-photometry recordings and optogenetic stimulation responses.
- **`04_Water_Consumption/`** — Analysis of behavioral events and event-aligned signals during water consumption.
- **`05_Open_Field_Test/`** — Open-field behavioral analysis, including locomotor patterns and movement-related measurements.
- **`06_Reaching_Test/`** — Quantification of reaching performance and DeepLabCut-based movement trajectories.

Each subfolder contains the corresponding scripts and a dedicated README with instructions and expected inputs and outputs.

## Data availability

The datasets associated with these analyses are available in the following repositories:

| Analysis | Data repository | Accession |
| --- | --- | --- |
| Bulk RNA-seq (`01_Bulk_RNAseq`) | NCBI GEO | **GSE233866**; additional bulk dataset: **GEOXXXXX** |
| MERFISH / MERSCOPE (`02_MERFISH`) | NCBI GEO | **GEOXXXXX** |
| Fiber photometry and optogenetics (`03_Fiber_Photometry_Optogenetics`) | DANDI Archive | [DANDI: 001828](https://doi.org/10.48324/dandi.001828/0.261020.1002) |
| Water consumption (`04_Water_Consumption`) | DANDI Archive | [DANDI: 001828](https://doi.org/10.48324/dandi.001828/0.261020.1002) |
| Open field test (`05_Open_Field_Test`) | DANDI Archive | [DANDI: 001828](https://doi.org/10.48324/dandi.001828/0.261020.1002) |
| Reaching test (`06_Reaching_Test`) | DANDI Archive | [DANDI: 001828](https://doi.org/10.48324/dandi.001828/0.261020.1002) |

## Usage

Refer to the README or configuration section within each subfolder. Set local input/output paths and install the dependencies listed in the relevant scripts before running the analyses.
