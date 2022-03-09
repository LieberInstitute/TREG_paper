library(SingleCellExperiment)
library(jaffelab)
library(edgeR)
library(limma)
library(tidyverse)
library(here)

#### Prep sce data ####
load(here("raw-data", "sce_pan.v2.Rdata"), verbose = TRUE)
summary(sce_pan$sum)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 102    5766    9813   17282   23563  196431 

#### Find Expression Rank ####
## takes like 50 min
rank_df <- apply(as.matrix(assays(sce_pan)$logcounts), 2, rank) %>%
  as.data.frame()

save(rank_df, file = here("processed-data", "01_find_tregs", "rank_df.Rdata"))

## subset to canididate genes "AKT3"   "ARID1B" "MALAT1" "POLR2A"
rank_df_subset <- rank_df[c("ENSG00000117020", "ENSG00000049618", "ENSG00000251562", "ENSG00000181222"),]
save(rank_df_subset, file = here("processed-data", "01_find_tregs", "rank_df_subset.Rdata"))

#### Linear modeling ####
mod_ct <- model.matrix(~log2(sum) + cellType.Broad, data = colData(sce_pan))

fit = lmFit(log2(assays(sce_pan)$counts + 1), mod_ct)
eB = eBayes(fit)
tt = topTable(eB, coef= "log2(sum)", number= Inf, sort.by = "none")
table(tt$adj.P.Val < 0.05)
# FALSE  TRUE 
#    14 23024 

save(tt, file = here("processed-data", "01_find_tregs", "lmfit.Rdata"))

# sgejobs::job_single('snRNA_validation', create_shell = TRUE, memory = '100G', command = "Rscript snRNA_validation.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

# [1] "Reproducibility information:"
# [1] "2022-03-09 14:32:25 EST"
# user    system   elapsed 
# 10514.456  2984.844 13510.703 
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-03-09
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# backports              1.4.1    2021-12-13 [2] CRAN (R 4.1.2)
# Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
# BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# broom                  0.7.12   2022-01-28 [2] CRAN (R 4.1.2)
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.1.0)
# cli                    3.2.0    2022-02-14 [2] CRAN (R 4.1.2)
# colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.1.2)
# crayon                 1.5.0    2022-02-14 [2] CRAN (R 4.1.2)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
# dbplyr                 2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# dplyr                * 1.0.8    2022-02-08 [2] CRAN (R 4.1.2)
# edgeR                * 3.36.0   2021-10-26 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  1.0.2    2022-01-14 [2] CRAN (R 4.1.2)
# forcats              * 0.5.1    2021-01-27 [2] CRAN (R 4.1.0)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
# gargle                 1.2.0    2021-07-02 [2] CRAN (R 4.1.0)
# generics               0.1.2    2022-01-31 [2] CRAN (R 4.1.2)
# GenomeInfoDb         * 1.30.1   2022-01-30 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
# GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
# ggplot2              * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
# glue                   1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
# googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# haven                  2.4.3    2021-08-04 [2] CRAN (R 4.1.1)
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# hms                    1.1.1    2021-09-26 [2] CRAN (R 4.1.2)
# httr                   1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# jaffelab             * 0.99.31  2021-10-29 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# jsonlite               1.8.0    2022-02-22 [2] CRAN (R 4.1.2)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [1] CRAN (R 4.1.1)
# limma                * 3.50.1   2022-02-17 [2] Bioconductor
# locfit                 1.5-9.5  2022-03-03 [2] CRAN (R 4.1.2)
# lubridate              1.8.0    2021-10-07 [2] CRAN (R 4.1.2)
# magrittr               2.0.2    2022-01-26 [2] CRAN (R 4.1.2)
# MASS                   7.3-55   2022-01-13 [3] CRAN (R 4.1.2)
# Matrix                 1.4-0    2021-12-08 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.1.0)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# pillar                 1.7.0    2022-02-01 [1] CRAN (R 4.1.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R6                     2.5.1    2021-08-19 [1] CRAN (R 4.1.1)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.1)
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.8    2022-01-13 [2] CRAN (R 4.1.2)
# RCurl                  1.98-1.6 2022-02-08 [2] CRAN (R 4.1.2)
# readr                * 2.1.2    2022-01-30 [2] CRAN (R 4.1.2)
# readxl                 1.3.1    2019-03-13 [2] CRAN (R 4.1.0)
# reprex                 2.0.1    2021-08-05 [2] CRAN (R 4.1.1)
# rlang                  1.0.2    2022-03-04 [2] CRAN (R 4.1.2)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.1.0)
# rvest                  1.0.2    2021-10-16 [2] CRAN (R 4.1.2)
# S4Vectors            * 0.32.3   2021-11-21 [2] Bioconductor
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# segmented              1.4-0    2022-01-28 [1] CRAN (R 4.1.2)
# sessioninfo            1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
# SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
# stringi                1.7.6    2021-11-29 [2] CRAN (R 4.1.2)
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
# tibble               * 3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
# tidyr                * 1.2.0    2022-02-01 [2] CRAN (R 4.1.2)
# tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
# tidyverse            * 1.3.1    2021-04-15 [2] CRAN (R 4.1.0)
# tzdb                   0.2.0    2021-10-27 [2] CRAN (R 4.1.2)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# withr                  2.5.0    2022-03-03 [2] CRAN (R 4.1.2)
# xml2                   1.3.3    2021-11-30 [2] CRAN (R 4.1.2)
# XVector                0.34.0   2021-10-26 [2] Bioconductor
# zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
# 
# [1] /users/lhuuki/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# **** Job ends ****
#   Wed Mar  9 14:36:22 EST 2022
# 
