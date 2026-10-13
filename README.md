# Anxa1_DAN_PD_codes

## Summary
This repository contains data and R scripts for analyzing the role of **Anxa1+ dopamine neurons (DANs) in Parkinson’s disease (PD)**. The experiments explore neural activity, behavioral effects of DAN silencing, and their impact on motor learning. Key findings include:
- Anxa1+ DANs are **selectively vulnerable** in early PD stages.
- Their activity correlates with **vigorous movements** rather than reward signaling.
- Silencing Anxa1+ DANs induces **bradykinesia-like motor deficits**.
- Anxa1+ DANs are essential for **procedural motor learning**.

## Description

This repository contains data, analysis code, and figures related to the experiments on Anxa1+ dopamine neurons (DANs) of the following preprint:

Mantas, I., Contestabile, A., Skara, V., Loiseau, C., Santos, I. A., Cramb, K. M. L., Filograna, R., Wade-Martins, R., Magill, P., & Meletis, K. (2024).  
Anxa1+ dopamine neuron vulnerability defines prodromal Parkinson’s disease bradykinesia and procedural motor learning impairment. bioRxiv.  
[https://doi.org/10.1101/2024.12.22.629963](https://doi.org/10.1101/2024.12.22.629963)

The datasets and scripts correspond to:

- **Figure 3**: Calcium and dopamine signal recordings during water consumption.
- **Figure 4**: Behavioral motif analysis and correlations with Anxa1+ DAN activity.
- **Figure 5**: Anxa1+ DAN silencing and its effect on motor behavior and motif usage.
- **Figure 6**: Procedural motor learning deficits following Anxa1+ DAN silencing.

## Repository Structure

```
/Anxa1_DAN_PD_experiments
│── /figure_3
│   ├── /data       # Raw and processed data for Figure 3D-E
│   ├── /code       # Scripts for analysis and visualization
│── /figure_4
│   ├── /data       # Raw and processed data for Figure 4C-H
│   ├── /code       # Scripts for analysis and visualization
│── /figure_5
│   ├── /data       # Raw and processed data for Figure 5B-M
│   ├── /code       # Scripts for analysis and visualization
│── /figure_6
│   ├── /data       # Raw and processed data for Figure 6B-D and 6H-J
│   ├── /code       # Scripts for analysis and visualization
│── README.md       # Repository documentation
│── LICENSE         # License file for usage
│── requirements.txt # Dependencies for analysis scripts
```

## Running Analysis
Each figure has a dedicated folder with **R** analysis scripts. To execute an analysis, use:
```bash
Rscript figure_3/code/analyze_fig3.R
```
Modify the script name accordingly for other figures.

## Data
The `data/` folders contain raw and processed datasets:
- **Raw data**: Includes behavioral, photometry, and optogenetics recordings.
- **Processed data**: Data matrices and extracted features used for analysis.

## Code
The `code/` folders contain all analysis scripts written in **R**.

## License
This repository is shared under the **CC-BY-NC-ND 4.0 License**, meaning it can be used for non-commercial purposes with proper attribution.
