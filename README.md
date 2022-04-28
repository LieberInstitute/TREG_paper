[![DOI](https://zenodo.org/badge/462890831.svg)](https://zenodo.org/badge/latestdoi/462890831)

# TREG paper

> Next-generation sequencing technologies have facilitated data-driven identification of gene sets with different features including genes with stable expression, cell-type specific expression, or spatially variable expression. Here, we aimed to define and identify a new class of "control" genes called Total RNA Expression Genes (TREGs), which correlate with total RNA abundance in heterogeneous cell types of different sizes and transcriptional activity. We provide a data-driven method to identify TREGs from single cell RNA-sequencing (RNA-seq) data, available as an R/Bioconductor package at https://bioconductor.org/packages/TREG. We demonstrated the utility of our method in the postmortem human brain using multiplex single molecule fluorescent in situ hybridization (smFISH) and compared candidate TREGs against classic housekeeping genes. We identified _AKT3_ as a top TREG across five brain regions, especially in the dorsolateral prefrontal cortex.

## Overview

In this repository we applied functions from `TREG` to find candidate TREGs in postmortem human brain snRNA-seq data from eight donors and five brain regions (Tran et al, Neuron, 2021 DOI: [10.1016/j.neuron.2021.09.001](https://doi.org/10.1016/j.neuron.2021.09.001)). We then validated top candidate TREGs (_AKT3_, _MALAT1_, and _ARID1B_) in different cell types of dorsolateral prefrontal cortex (DLPFC) using smFISH with RNAscope technology then analyzed with HALO.

## Cite this work

If you use anything in this repository please cite the following publication:

Huuki-Myers LA, Montgomery KD, Kwon SH, Page SC, Hicks SC, Maynard KR, Collado-Torres L (2022). “Data Driven Identification of Total RNA Expression Genes "TREGs" for estimation of RNA abundance in heterogeneous cell types.” _bioRxiv_. doi: [10.1101/TODO](https://doi.org/10.1101/TODO).

A BibTeX entry for LaTeX users is

```
@Article{TREGpaper,
    title = {Data-driven Identification of Total RNA Expression Genes (TREGs) for Estimation of RNA Abundance in Heterogeneous Cell Types},
    author = {Louise A. Huuki-Myers and Kelsey D. Montgomery and Sang Ho. Kwon and Stephanie C. Page and Stephanie C. Hicks and Kristen R. Maynard and Leonardo Collado-Torres},
    year = {2022},
    journal = {bioRxiv},
    doi = {10.1101/TODO},
    url = {https://doi.org/10.1101/TODO},
  }
```

## Organization

Organization of this project is guided by [R/Bioconductor-powered Team Data Science](https://lcolladotor.github.io/bioc_team_ds/organizing-your-work.html#.YkNW_DfMJfV.)

* `raw-data`: HALO output, RNAscope images, note most of these files are too big for GitHub but can be accessed via [Globus](http://research.libd.org/globus/)

* `code`: R and shell scripts
  + `01_find_tregs`: Run the TREG discovery pipeline, explore results in the snRNA-seq data  
  + `02_analyze_halo`: Import, quality control, and analyze halo segmentation data  from smFISH validation 

* `processed-data`: Data output from analysis as `.Rdata` or `.csv` files

* `plots`: pdfs used in main & supplementary figures are separated from other exploratory plots created in this project. 

## Data availability

The raw data for this project and all the results are publicly available through Globus at [jhpce#TREG_paper](http://research.libd.org/globus/jhpce_TREG_paper/index.html). 

While every computer system is different, you might benefit from checking the [using Globus to transfer files](https://jhpce.jhu.edu/knowledge-base/using-globus-to-transfer-files/) tutorial for [JHPCE](https://jhpce.jhu.edu/) users.

## License

<img src="https://licensebuttons.net/l/by-nc/3.0/88x31.png" alt width="88" height="31" scale="0">
Attribution-NonCommercial: CC BY-NC

This license lets others remix, tweak, and build upon our work non-commercially as long as they acknowledge our work.

[View License Deed](https://creativecommons.org/licenses/by-nc/4.0) | [View Legal Code](https://creativecommons.org/licenses/by-nc/4.0/legalcode)

## Internal

JHPCE location: `/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/TREG_paper`.
