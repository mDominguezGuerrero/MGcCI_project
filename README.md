Bioinformatic analysis of “single-nuclei RNA sequencing” from mouse
brain with microglia deficient in mitochondrial complex I
================
Máximo Domínguez Guerrero, Alberto Pascual Bravo, Juan José Pérez Moreno

This repository contains the code and documentation for the
bioinformatic analysis of “single-nuclei RNA sequencing” from mouse
brain with microglia deficient in mitochondrial complex I.

Data can be accessed on request to the authors.

# Repository Structure

- **R/**: This folder contains functions developed for analysis.
- **MGcCI.md** and **MGcCI.Rmd**: The main markdown file with all
  analyses.

# Getting Started

## Prerequisites

- R (version 4.3.2 or later)
- R packages:
- AnnotationDbi
- ComplexHeatmap
- EnhancedVolcano
- HGNChelper
- Seurat
- VennDetail
- dplyr
- garnett
- ggplot2
- ggpubr
- ggVennDiagram
- harmony
- monocle3
- openxlsx
- org.Mm.eg.db

## Usage

1.  Load the functions from the **R** folder:

``` r
source("R/assign_group.R")
source("R/extract_module.R")
source("R/filter_cells.R")
```

2.  Open the **MCcCI.md** or **MGcCI.Rmd** file to follow along with the
    analyses.