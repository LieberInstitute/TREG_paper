[![DOI](https://zenodo.org/badge/462890831.svg)](https://zenodo.org/badge/latestdoi/462890831)

# TREG paper

## Overview
We aimed to define and identify a new class of "control" genes called Total RNA Expression Genes (TREGs), which correlate with total RNA abundance in heterogeneous cell types of different sizes and transcriptional activity. We developed softeware to identify TREGs from single nucleus RNA-seq data, available as an R/Bioconductor `TREG` package at http://research.libd.org/TREG/.

In this repository we applied functions from `TREG` to find candidate TREGs in postmortem human brain snRNA-seq data from eight donors and five brain regions (Tran et al, Neuron, 2021). We then validated top TREGs (AKT3, MALAT1, and ARID1B) in different cell types of dorsolateral prefrontal cortex (DLPFC) using smFISH with RNAscope technology then analyzed with HALO.

### Organization
Organization of this project is guided by [R/Bioconductor-powered Team Data Science](https://lcolladotor.github.io/bioc_team_ds/organizing-your-work.html#.YkNW_DfMJfV.)

* `raw-data`: HALO output, RNAscope images, note most of these files are too big for github but can be accessed via [globus](http://research.libd.org/globus/)

* `code`: R and shell scripts
  + `01_find_tregs`: Run the TREG discovery pipeline, explore results in the snRNA-seq data  
  + `02_analyze_halo`: Import, quality control, and analyze halo segmentation data  from smFISH validation 

* `processed-data`: Data output from analysis as `.Rdata` or `.csv` files

* `plots`: pdfs used in main & supplementary figures are separated from other exploratory plots created in this project. 

## Cite this work

TODO

## Internal

JHPCE location: `/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/TREG_paper`.
