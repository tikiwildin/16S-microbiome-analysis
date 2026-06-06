# # 16S rRNA Microbiome Analysis
This repository contains R scripts and example figures from a 16S rRNA microbiome analysis project focused on soil- and plant-associated microbial communities. The workflow includes taxonomic visualization, alpha diversity analysis, beta diversity analysis, and predicted functional profiling using PICRUSt2 outputs.
## Project Overview

The goal of this project was to analyze microbial community composition across different sample types and experimental conditions, including native ant mound soil, cultured ant mound samples, and wheatgrass-associated microbiomes.

The analyses in this repository focus on:

- Genus-level taxonomic composition
- Shannon alpha diversity
- Bray-Curtis beta diversity and PCoA visualization
- Predicted functional pathway profiles using PICRUSt2

## Repository Structure

```text
16S-microbiome-analysis/
│
├── figures/
│   ├── antmound_genus_stacked_barplot.png
│   ├── BetaDiversityPCoA.png
│   ├── picrust2_wheatgrass_heatmap_chronological.png
│   └── ShannonBoxplot.png
│
├── scripts/
│   ├── AlphaDiversityShannonboxplot.R
│   ├── BetaDiversityBrayCurtisForAntmound.R
│   ├── picrustpathwayanalysisforWheatgrass.R
│   └── plot_antmound_genus_barplot.R
│
└── README.md
```
## Scripts

AlphaDiversityShannonboxplot.R:
- Generates Shannon alpha diversity boxplots for comparing microbial diversity across sample groups.
<img width="600" height="400" alt="image" src="https://github.com/user-attachments/assets/68d4681d-aa76-46b7-83f2-fe2ab89208fe" />

BetaDiversityBrayCurtisForAntmound.R:
- Reads a mothur generated Bray-Curtis distance matrix, performs Principle Coordinates Analysis, and visulizes beta diversity patterns among ant mound samples
<img width="700" height="400" alt="BetaDiversityPCoA" src="https://github.com/user-attachments/assets/45748c68-6da9-4aa2-bf08-8913a9a9f18f" />

picrustpathwayanalysisforWheatgrass.R
- Analyzes PICRUSt2 predicted pathway abundance data for wheatgrass samples. The script normalizes pathway abundances, performs PCA, and generates a chronological heatmap for the most variable predicted pathways.
<img width="600" height="760" alt="picrust2_wheatgrass_heatmap_chronological" src="https://github.com/user-attachments/assets/c2e7fd43-e500-41e2-a232-6fdb75ed8db3" />

plot_antmound_genus_barplot.R
- Reads mothur/Krona genus-level taxonomy files, calculates relative abundance, and generates a stacked barplot showing dominant genera across ant mound sample groups
<img width="750" height="400" alt="antmound_genus_stacked_barplot" src="https://github.com/user-attachments/assets/bfd2e1aa-b3a3-48a4-b638-a5a7152cf9f0" />

## Tools and Packages Used
- R
- tidyverse
- ggplot2
- dplyr
- ape
- pheatmap
- PICRUSt2
- mothur/Galaxy-based 16S processing outputs

## Data Availability
Raw sequencing data and lab-specific metadata are not included in this repository because they are part of an unpublished research project. The scripts are provided to demonstrate the analysis workflow and organization of the computational methods.

## Author
Atiq Bacus
Hunter College
