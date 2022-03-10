library(SingleCellExperiment)
library(SummarizedExperiment)
library(TREG)
library(tidyverse)
library(here)
library(sessioninfo)

load(here("raw-data", "sce_pan.v2.Rdata"), verbose = TRUE)

#### Refine Region Lables ####
## Capitalize Regions
sce_pan$region2 <- toupper(sce_pan$region)
sce_pan$region2[sce_pan$region2 == "SACC"] <- "sACC"
sce_pan$region2[sce_pan$region2 == "NAC"] <- "NAc"

## Examine the count for each cell type
(nuclei_counts <- table(sce_pan$cellType.Broad, sce_pan$region2))
#         AMY DLPFC   HPC   NAc  sACC
# Astro  1638   782  1170  1099   907
# Endo     31     0     0     0     0
# Macro     0    10     0    22     0
# Micro  1168   388  1126   492   784
# Mural    39    18    43     0     0
# Oligo  6080  5455  5912  6134  4584
# OPC    1459   572   838   669   911
# Tcell    31     9    26     0     0
# Excit   443  2388   623     0  4163
# Inhib  3117  1580   366 11476  3974

write.csv(nuclei_counts, file = here("processed-data", "01_find_tregs", "supp_tables", "nuclei_counts.csv"))

## Rare cell types
sce_pan$region2[sce_pan$cellType.Broad %in% c("Macro", "Mural", "Endo", "Tcell")] <- "combined"
table(sce_pan$region2)

sce_pan$ctXregion <- paste0(sce_pan$cellType.Broad, "_", sce_pan$region2)

## Calc sum_counts
# sce_pan$sum_counts <- colSums(assays(sce_pan)$counts)

#### filter for top 50% expressed of genes ####
## record data for each gene
gene_metrics <- rowData(sce_pan) %>%
    as.data.frame() %>%
    select(Symbol = gene_name, ensembl_id = gene_id)

dim(sce_pan)
# [1] 23038 70527
row_means <- rowMeans(assays(sce_pan)$logcounts)
median_row_means <- median(row_means)

sce_pan <- sce_pan[row_means > median_row_means, ]
dim(sce_pan)
# [1] 11519 70527

gene_metrics$top50 <- gene_metrics$ensembl_id %in% rownames(sce_pan)

#### Proportion Zero Filtering ####
gene_propZero <- get_prop_zero(sce_pan, "ctXregion")
head(gene_propZero)

## record max prop zero for gene metrics
max_propZero <- apply(gene_propZero, 1, max)
gene_metrics$max_PropZero <- NA
gene_metrics[names(max_propZero), ]$max_PropZero <- max_propZero

## filter by Proportion Zero
propZero_limit <- 0.75
genes_filtered <- filter_prop_zero(gene_propZero, cutoff = propZero_limit)

# max_prop_zero <- apply(gene_propZero, 1, max)
# genes_filtered <- names(max_prop_zero[max_prop_zero < propZero_limit])
length(genes_filtered)
# [1] 877

sce_pan <- sce_pan[genes_filtered, ]
dim(sce_pan)
# [1]   877 70527

gene_metrics$PropZero_filter <- gene_metrics$ensembl_id %in% rownames(sce_pan)

#### Run Rank Invariance for full data set ####
assays(sce_pan)$logcounts <- as.matrix(assays(sce_pan)$logcounts)

rank_invariance <- rank_invariance_express(sce_pan, group_col = "cellType.Broad")

rank_invar_df <- as.data.frame(rank_invariance)
colnames(rank_invar_df) <- "rank_invar"

## add to gene metrics
gene_metrics <- gene_metrics %>%
    left_join(rank_invar_df %>%
        rownames_to_column("ensembl_id"))

head(gene_metrics)

gene_metrics %>%
    arrange(-rank_invar) %>%
    head(10)

#        Symbol      ensembl_id top50 max_PropZero PropZero_filter rank_invar
# 1      MALAT1 ENSG00000251562  TRUE  0.004065041            TRUE        877
# 2      JMJD1C ENSG00000171988  TRUE  0.212121212            TRUE        876
# 3         FTX ENSG00000230590  TRUE  0.212121212            TRUE        875
# 4       MACF1 ENSG00000127603  TRUE  0.156250000            TRUE        874
# 5        AKT3 ENSG00000117020  TRUE  0.531250000            TRUE        873
# 6      TNRC6B ENSG00000100354  TRUE  0.218750000            TRUE        872
# 7  AC016831.7 ENSG00000285106  TRUE  0.225806452            TRUE        871
# 8      ZFAND3 ENSG00000156639  TRUE  0.198734729            TRUE        870
# 9      ARID1B ENSG00000049618  TRUE  0.322580645            TRUE        869
# 10      CADM2 ENSG00000175161  TRUE  0.734458259            TRUE        868

## Write gene_metrics to csv for Supplement
write.csv(gene_metrics, file = here("processed-data", "01_find_tregs", "supp_tables", "gene_metrics.csv"))

## Extract RI data for analysis
rank_invar_df <- gene_metrics %>%
    filter(!is.na(rank_invar)) %>%
    select(Symbol, ensembl_id, rank_invar) %>%
    arrange(-rank_invar)

head(rank_invar_df, 10)
# Symbol      ensembl_id rank_invar
# 1      MALAT1 ENSG00000251562        877
# 2      JMJD1C ENSG00000171988        876
# 3         FTX ENSG00000230590        875
# 4       MACF1 ENSG00000127603        874
# 5        AKT3 ENSG00000117020        873
# 6      TNRC6B ENSG00000100354        872
# 7  AC016831.7 ENSG00000285106        871
# 8      ZFAND3 ENSG00000156639        870
# 9      ARID1B ENSG00000049618        869
# 10      CADM2 ENSG00000175161        868

save(gene_propZero, rank_invar_df, file = here("processed-data", "01_find_tregs", "supp_tables", "rank_invar.Rdata"))

# sgejobs::job_single('find_TREGs', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript find_TREGs.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-02-24
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
#
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# backports              1.4.1    2021-12-13 [2] CRAN (R 4.1.2)
# Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
# BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# broom                  0.7.12   2022-01-28 [2] CRAN (R 4.1.2)
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.1.0)
# cli                    3.2.0    2022-02-14 [2] CRAN (R 4.1.2)
# codetools              0.2-18   2020-11-04 [3] CRAN (R 4.1.2)
# colorout             * 1.2-2    2021-12-03 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.1.2)
# crayon                 1.5.0    2022-02-14 [2] CRAN (R 4.1.2)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
# dbplyr                 2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# dplyr                * 1.0.8    2022-02-08 [2] CRAN (R 4.1.2)
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  1.0.2    2022-01-14 [2] CRAN (R 4.1.2)
# forcats              * 0.5.1    2021-01-27 [2] CRAN (R 4.1.0)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
# generics               0.1.2    2022-01-31 [2] CRAN (R 4.1.2)
# GenomeInfoDb         * 1.30.1   2022-01-30 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
# GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
# ggplot2              * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
# glue                   1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# haven                  2.4.3    2021-08-04 [2] CRAN (R 4.1.1)
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# hms                    1.1.1    2021-09-26 [2] CRAN (R 4.1.2)
# httr                   1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# jsonlite               1.8.0    2022-02-22 [2] CRAN (R 4.1.2)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [1] CRAN (R 4.1.1)
# lubridate              1.8.0    2021-10-07 [2] CRAN (R 4.1.2)
# magrittr               2.0.2    2022-01-26 [2] CRAN (R 4.1.2)
# Matrix                 1.4-0    2021-12-08 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.1.0)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# pillar                 1.7.0    2022-02-01 [1] CRAN (R 4.1.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# pryr                   0.1.5    2021-07-26 [2] CRAN (R 4.1.0)
# purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R6                     2.5.1    2021-08-19 [1] CRAN (R 4.1.1)
# rafalib                1.0.0    2015-08-09 [1] CRAN (R 4.1.1)
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.8    2022-01-13 [2] CRAN (R 4.1.2)
# RCurl                  1.98-1.6 2022-02-08 [2] CRAN (R 4.1.2)
# readr                * 2.1.2    2022-01-30 [2] CRAN (R 4.1.2)
# readxl                 1.3.1    2019-03-13 [2] CRAN (R 4.1.0)
# reprex                 2.0.1    2021-08-05 [2] CRAN (R 4.1.1)
# rlang                  1.0.1    2022-02-03 [2] CRAN (R 4.1.2)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.1.0)
# rvest                  1.0.2    2021-10-16 [2] CRAN (R 4.1.2)
# S4Vectors            * 0.32.3   2021-11-21 [2] Bioconductor
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
# sgejobs                0.99.1   2021-12-03 [1] Github (LieberInstitute/sgejobs@f5ab0ca)
# SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
# stringi                1.7.6    2021-11-29 [2] CRAN (R 4.1.2)
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
# tibble               * 3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
# tidyr                * 1.2.0    2022-02-01 [2] CRAN (R 4.1.2)
# tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
# tidyverse            * 1.3.1    2021-04-15 [2] CRAN (R 4.1.0)
# TREG                 * 0.99.0   2022-02-24 [1] Github (LieberInstitute/TREG@d6106d6)
# tzdb                   0.2.0    2021-10-27 [2] CRAN (R 4.1.2)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# withr                  2.4.3    2021-11-30 [2] CRAN (R 4.1.2)
# xml2                   1.3.3    2021-11-30 [2] CRAN (R 4.1.2)
# XVector                0.34.0   2021-10-26 [2] Bioconductor
# zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
#
# [1] /users/lhuuki/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
#
# ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#
