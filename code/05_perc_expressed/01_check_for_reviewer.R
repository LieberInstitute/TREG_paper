library("SingleCellExperiment")
library("ggplot2")
library("here")
library("sessioninfo")

## Plot output directories
dir_plots <- here::here(
    "plots",
    "05_perc_expressed"
)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

## Load the snRNA-seq data
load(here("raw-data", "sce_pan.v2.Rdata"), verbose = TRUE)

dim(sce_pan)
# [1] 23038 70527
nuclei_expressed <- rowSums(counts(sce_pan) > 0)
summary(nuclei_expressed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  100    1191    5858   12464   19662   70522

perc_expressed <- nuclei_expressed / ncol(sce_pan) * 100
summary(perc_expressed)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1418  1.6891  8.3053 17.6733 27.8780 99.9929

expr_above_perc <- sapply(seq(0, 100, by = 1), function(x) {
    sum(perc_expressed >= x)
})

## Put it together in a dataframe
df <- data.frame(perc = seq(0, 100, by = 1), expr_above_perc)
df
#     perc expr_above_perc
# 1      0           23038
# 2      1           18997
# 3      2           16726
# 4      3           15383
# 5      4           14364
# 6      5           13541
# 7      6           12829
# 8      7           12217
# 9      8           11682
# 10     9           11207
# 11    10           10771
# 12    11           10369
# 13    12           10003
# 14    13            9634
# 15    14            9314
# 16    15            9002
# 17    16            8696
# 18    17            8427
# 19    18            8143
# 20    19            7866
# 21    20            7603
# 22    21            7323
# 23    22            7077
# 24    23            6846
# 25    24            6606
# 26    25            6378
# 27    26            6138
# 28    27            5933
# 29    28            5724
# 30    29            5511
# 31    30            5303
# 32    31            5091
# 33    32            4900
# 34    33            4720
# 35    34            4542
# 36    35            4355
# 37    36            4214
# 38    37            4038
# 39    38            3871
# 40    39            3722
# 41    40            3592
# 42    41            3419
# 43    42            3265
# 44    43            3135
# 45    44            3005
# 46    45            2879
# 47    46            2745
# 48    47            2607
# 49    48            2500
# 50    49            2381
# 51    50            2255
# 52    51            2157
# 53    52            2046
# 54    53            1940
# 55    54            1825
# 56    55            1736
# 57    56            1667
# 58    57            1587
# 59    58            1514
# 60    59            1441
# 61    60            1369
# 62    61            1301
# 63    62            1228
# 64    63            1171
# 65    64            1099
# 66    65            1032
# 67    66             981
# 68    67             913
# 69    68             868
# 70    69             809
# 71    70             769
# 72    71             727
# 73    72             669
# 74    73             620
# 75    74             567
# 76    75             527
# 77    76             487
# 78    77             446
# 79    78             410
# 80    79             370
# 81    80             338
# 82    81             299
# 83    82             260
# 84    83             243
# 85    84             220
# 86    85             196
# 87    86             172
# 88    87             151
# 89    88             121
# 90    89             105
# 91    90              92
# 92    91              75
# 93    92              61
# 94    93              45
# 95    94              31
# 96    95              22
# 97    96              15
# 98    97               6
# 99    98               2
# 100   99               1
# 101  100               0

pdf(file.path(dir_plots, "genes_percent_expressed.pdf"), height = 8, width = 8)
ggplot(
    df,
    aes(x = perc, y = expr_above_perc)
) +
    geom_point() +
    geom_line() +
    geom_smooth() +
    theme_bw(base_size = 20) +
    ylab("N genes expressed >= percent of nuclei") +
    xlab("Percent of nuclei where this gene is expressed")
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.1 Patched (2023-06-16 r84551)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2023-06-16
#  pandoc   3.1.1 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
#  BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
#  bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.0)
#  cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.0)
#  colorout               1.2-2     2023-05-06 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.0)
#  crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.0)
#  DelayedArray           0.26.3    2023-05-22 [2] Bioconductor
#  digest                 0.6.31    2022-12-11 [2] CRAN (R 4.3.0)
#  dplyr                  1.1.2     2023-04-20 [2] CRAN (R 4.3.0)
#  fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.0)
#  farver                 2.1.1     2022-07-06 [2] CRAN (R 4.3.0)
#  fastmap                1.1.1     2023-02-24 [2] CRAN (R 4.3.0)
#  generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.0)
#  GenomeInfoDb         * 1.36.0    2023-04-25 [2] Bioconductor
#  GenomeInfoDbData       1.2.10    2023-04-11 [2] Bioconductor
#  GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
#  ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.3.0)
#  glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.0)
#  gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.3.0)
#  here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.0)
#  htmltools              0.5.5     2023-03-23 [2] CRAN (R 4.3.0)
#  htmlwidgets            1.6.2     2023-03-17 [2] CRAN (R 4.3.0)
#  httpuv                 1.6.11    2023-05-11 [2] CRAN (R 4.3.0)
#  IRanges              * 2.34.0    2023-04-25 [2] Bioconductor
#  jsonlite               1.8.5     2023-06-05 [2] CRAN (R 4.3.0)
#  labeling               0.4.2     2020-10-20 [2] CRAN (R 4.3.0)
#  later                  1.3.1     2023-05-02 [2] CRAN (R 4.3.0)
#  lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
#  lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.0)
#  magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.0)
#  Matrix                 1.5-4.1   2023-05-18 [3] CRAN (R 4.3.1)
#  MatrixGenerics       * 1.12.2    2023-06-09 [2] Bioconductor
#  matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.0)
#  mgcv                   1.8-42    2023-03-02 [3] CRAN (R 4.3.1)
#  mime                   0.12      2021-09-28 [2] CRAN (R 4.3.0)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.0)
#  nlme                   3.1-162   2023-01-31 [3] CRAN (R 4.3.1)
#  pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.0)
#  pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.0)
#  png                    0.1-8     2022-11-29 [2] CRAN (R 4.3.0)
#  promises               1.2.0.1   2021-02-11 [2] CRAN (R 4.3.0)
#  R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.0)
#  Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.3.0)
#  RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.0)
#  rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.0)
#  rmote                  0.3.4     2023-05-06 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.0)
#  S4Arrays               1.0.4     2023-05-14 [2] Bioconductor
#  S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
#  scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.0)
#  servr                  0.27      2023-05-02 [1] CRAN (R 4.3.0)
#  sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.0)
#  SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
#  SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
#  tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.0)
#  tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.0)
#  utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.0)
#  vctrs                  0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
#  withr                  2.5.0     2022-03-03 [2] CRAN (R 4.3.0)
#  xfun                   0.39      2023-04-20 [2] CRAN (R 4.3.0)
#  XVector                0.40.0    2023-04-25 [2] Bioconductor
#  zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
# 
#  [1] /users/lcollado/R/4.3
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/R/4.3/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/R/4.3/lib64/R/library
# 
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
