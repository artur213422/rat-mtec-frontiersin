# Conserved and species-specific Aire-dependent programs in thymic epithelial cells

This repository contains the analysis code and processed inputs used to generate figure panels for the manuscript:

**Conserved and species-specific Aire-dependent programs in thymic epithelial cells**

The workflow is organized in Python, with R scripts used for statistical analyses. 
Python handles data preparation, processing, and visualization, primarily using Polars and Matplotlib. 
Running the upstream analysis workflow requires both Python and R.

The code in `src/manuscript/` generates manuscript figure panels from the included processed data. 
Panels are exported primarily as SVG files and were subsequently edited manually for publication. 
The generated panels therefore differ in appearance from the final published figures.

The remaining code documents the upstream processing and analysis workflow. 
Running `python bootstrap.py` creates the expected directory structure and downloads a subset of the required source data 
and reference resources. Additional inputs must be obtained separately, including the CEL files described in the manuscript. 


## Repository structure

```text
config/                  Dataset definitions and processing settings
data/                    Selected standardized data and manuscript inputs
src/
├── manuscript/          Manuscript figure and panel generation
├── analysis/            Upstream analyses
├── data_preprocessing/  Source-data processing
├── resource_processing/ Reference-data preparation
└── common/              Shared functions and utilities
bootstrap.py             Directory setup and available input downloads
```

The main entry points for figure generation are in `src/manuscript/figures/`. Individual panel scripts are in `src/manuscript/panels/`.
