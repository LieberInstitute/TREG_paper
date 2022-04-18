library("SingleCellExperiment")
library("TREG")
library("tidyverse")
library("jaffelab")
library("here")
library("sessioninfo")
library("GGally")
library("UpSetR")

## Load data
plot_dir <- "plots/01_find_tregs"
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)
load(here("raw-data", "sce_pan.v2.Rdata"), verbose = TRUE)
gene_metrics <- read.csv(here("processed-data", "01_find_tregs", "supp_tables", "gene_metrics.csv"), row.names = 1)

#### Prep Data ####
## Refine region lables
sce_pan$region <- toupper(sce_pan$region)
sce_pan$region[sce_pan$region == "SACC"] <- "sACC"
sce_pan$region[sce_pan$region == "NAC"] <- "NAc"

table(sce_pan$region)
table(sce_pan$region, sce_pan$cellType.Broad)
#   AMY DLPFC   HPC   NAc  sACC
# 14006 11202 10104 19892 15323

sce_pan$region2 <- sce_pan$region
sce_pan$region2[sce_pan$cellType.Broad %in% c("Macro", "Mural", "Endo", "Tcell")] <- "combined"
sce_pan$region <- as.factor(sce_pan$region)
sce_pan$region2 <- as.factor(sce_pan$region2)
#
# sce_pan$ctXregion <- paste0(sce_pan$cellType.Broad, "_", sce_pan$region2)
# table(sce_pan$ctXregion)

#### Run Full filter + RI process per region ####
## filter to top 50% expression & run get_prop_zero()
region_propZero <- map2(splitit(sce_pan$region), levels(sce_pan$region), function(region_i, region_name) {
    Sys.time()
    sce_region <- sce_pan[, region_i]
    message("Region: ", region_name, ", n_nuc = ", ncol(sce_region))
    row_means <- rowMeans(assays(sce_region)$logcounts)
    median_row_means <- median(row_means)
    sce_region <- sce_region[row_means > median_row_means, ]
    message("Median Express: ", round(median_row_means, 3), ", n_gene = ", nrow(sce_region))
    ## get prop zero
    gene_propZero <- get_prop_zero(sce_region, "cellType.Broad")
    return(gene_propZero)
})

# Region: AMY, n_nuc = 14006
# Median Express: 0.047, n_gene = 11519
# Region: DLPFC, n_nuc = 11202
# Median Express: 0.053, n_gene = 11519
# Region: HPC, n_nuc = 10104
# Median Express: 0.047, n_gene = 11519
# Region: NAc, n_nuc = 19892
# Median Express: 0.059, n_gene = 11519
# Region: sACC, n_nuc = 15323
# Median Express: 0.053, n_gene = 11519

## Plot Prop Zero densities
pdf(here(plot_dir, "supp_pdf", "region_propZero_density.pdf"), width = 10)
map2(region_propZero, names(region_propZero), function(propZero, region_name) {
    propZero_long <- propZero %>%
        rownames_to_column("gene") %>%
        pivot_longer(!gene, values_to = "propZero") %>%
        filter(!is.na(propZero)) %>%
        mutate(cellType = factor(name, levels = levels(sce_pan$cellType.Broad)))

    propZero_long %>%
        ggplot(aes(x = propZero, fill = cellType)) +
        geom_histogram(binwidth = 0.05, color = "black", size = .2) +
        scale_fill_manual(values = cell_colors) +
        facet_wrap(~cellType) +
        labs(
            title = region_name,
            x = "Proportion Zero per Group",
            y = "Number of Genes"
        ) +
        geom_vline(xintercept = 0.75, color = "red", linetype = "dashed") +
        theme_bw() +
        theme(
            legend.position = "none",
            text = element_text(size = 15),
            axis.text.x = element_text(angle = 90, hjust = 1)
        ) +
        scale_y_continuous(breaks = seq(0, 3000, 1500)) +
        scale_x_continuous(breaks = seq(0, 1, .5))
})
dev.off()

#### Calc RI build gene metrics for each region ####

gene_metrics_region <- map2(region_propZero, splitit(sce_pan$region), function(propZero, region_i) {
    ## Annotate top 50% genes
    gene_metrics_r <- gene_metrics_blank
    gene_metrics_r$top50 <- gene_metrics_r$ensembl_id %in% rownames(propZero)

    ## filter propZero
    ## record max prop zero for gene metrics
    max_propZero <- apply(propZero, 1, max, na.rm = TRUE)
    gene_metrics_r$max_PropZero <- NA
    gene_metrics_r[names(max_propZero), ]$max_PropZero <- max_propZero

    ## filter by Proportion Zero
    propZero_limit <- 0.75
    genes_filtered <- filter_prop_zero(propZero, cutoff = propZero_limit)
    # return(genes_filtered)
    message("filter to: ", length(genes_filtered))

    sce_pan_temp <- sce_pan[genes_filtered, ]
    sce_pan_temp <- sce_pan_temp[, region_i]
    message(str(dim(sce_pan_temp)))

    gene_metrics_r$PropZero_filter <- gene_metrics_r$ensembl_id %in% genes_filtered

    #### Run Rank Invariance for full data set ####
    assays(sce_pan_temp)$logcounts <- as.matrix(assays(sce_pan_temp)$logcounts)

    rank_invariance <- rank_invariance_express(sce_pan_temp, group_col = "cellType.Broad")

    rank_invar_df <- as.data.frame(rank_invariance)
    colnames(rank_invar_df) <- "rank_invar"

    ## add to gene metrics
    gene_metrics_r <- gene_metrics_r %>%
      left_join(rank_invar_df %>%
                  rownames_to_column("ensembl_id"),  by = "ensembl_id")

    return(gene_metrics_r)
})

gene_metrics_all <- c(list(ALL = gene_metrics), gene_metrics_region)

map(gene_metrics_region, ~.x %>% arrange(-rank_invar) %>% head())
map(gene_metrics_region, ~.x %>% mutate(r = sum(!is.na(.x$rank_invar)) - rank_invar + 1) %>%filter(Symbol %in% c("MALAT1", "AKT3", "ARID1B")))

## get top 10% of RI genes, make upset plots
top_ri <- map(gene_metrics_region, ~.x %>% 
                     arrange(-rank_invar) %>%
                     # head(.x %>% 
                     #          filter(!is.na(rank_invar)) %>%
                     #          nrow()%/%10) %>%
                     head(40) %>%
                     pull(Symbol))


pdf(here(plot_dir, "supp_pdf", "upset_region_top_ri.pdf"))
upset(fromList(top_ri), order.by = "freq")
dev.off()

Reduce(intersect, top_ri)
# [1] "MALAT1"     "MACF1"      "AKT3"       "TNRC6B"     "JMJD1C"     "FTX"        "AC016831.7" "ZFAND3"    
# [9] "ARID1B"     "KMT2C"      "RERE"       "KANSL1"     "MED13L" 

head(gene_metrics_region[[1]])
## filter to 877 Prop Zero genes
gene_set <- gene_metrics$ensembl_id[gene_metrics$PropZero_filter]
length(gene_set)
# [1] 877

sce_pan <- sce_pan[gene_set, ]
dim(sce_pan)
# [1]   877 70527

#### Run Rank Invariance for full data set ####

region_RI <- map_dfc(
    splitit(sce_pan$region),
    ~ rank_invariance_express(sce_pan[, .x], group_col = "cellType.Broad")
)

ctXregion_RI <- rank_invariance_express(sce_pan, group_col = "ctXregion")

## join w/ gene metrics
gene_metrics2 <- gene_metrics %>%
    left_join(region_RI %>% add_column(
        ensembl_id = gene_set,
        ctXregion = ctXregion_RI
    ))

gene_metrics2 %>%
    arrange(-rank_invar) %>%
    head(10)
#        Symbol      ensembl_id top50 max_PropZero PropZero_filter rank_invar AMY DLPFC HPC NAc  sACC
# 1      MALAT1 ENSG00000251562  TRUE  0.004065041            TRUE        877 877   877 877 877 877.0
# 5        AKT3 ENSG00000117020  TRUE  0.531250000            TRUE        873 875   869 874 859 869.0
# 9      ARID1B ENSG00000049618  TRUE  0.322580645            TRUE        869 869   870 876 868 872.0

## Write gene_metrics to csv for Supplement
write.csv(gene_metrics2, file = here("processed-data", "01_find_tregs", "supp_tables", "gene_metrics_regions.csv"))


#### Plots ####
plot_dir <- "plots/01_find_tregs"

scatter_RI <- ggpairs(gene_metrics2 %>% filter(PropZero_filter),
    columns = c("rank_invar", "AMY", "DLPFC", "HPC", "NAc", "sACC", "ctXregion"),
    lower = list(continuous = wrap("points", alpha = 0.4, size = 0.15))
)

ggsave(scatter_RI, filename = here(plot_dir, "explore", "region_RI_scatter.png"), height = 10, width = 10)


scatter_RI_top10p <- ggpairs(gene_metrics2 %>% filter(rank_invar > 850),
    columns = c("rank_invar", "AMY", "DLPFC", "HPC", "NAc", "sACC", "ctXregion")
    # ,
    # lower = list(continuous = wrap("points", alpha = 0.4,size=0.15))
)

ggsave(scatter_RI_top10p, filename = here(plot_dir, "explore", "region_RI_top10p_scatter.png"), height = 10, width = 10)



# sgejobs::job_single('find_TREGs', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript find_TREGs.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R Under development (unstable) (2021-11-06 r81149)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-04-06
# pandoc   2.11.0.4 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/bin/pandoc
#
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# backports              1.4.1    2021-12-13 [2] CRAN (R 4.2.0)
# Biobase              * 2.55.1   2022-03-31 [2] Bioconductor
# BiocGenerics         * 0.41.2   2021-11-15 [2] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.2.0)
# broom                  0.7.12   2022-01-28 [2] CRAN (R 4.2.0)
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.1.0)
# cli                    3.2.0    2022-02-14 [2] CRAN (R 4.2.0)
# colorout             * 1.2-2    2022-04-06 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.2.0)
# crayon                 1.5.1    2022-03-26 [2] CRAN (R 4.2.0)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.2.0)
# dbplyr                 2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
# DelayedArray           0.21.2   2021-11-16 [2] Bioconductor
# dplyr                * 1.0.8    2022-02-08 [2] CRAN (R 4.2.0)
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.2.0)
# fansi                  1.0.3    2022-03-24 [2] CRAN (R 4.2.0)
# forcats              * 0.5.1    2021-01-27 [2] CRAN (R 4.1.0)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.2.0)
# gargle                 1.2.0    2021-07-02 [2] CRAN (R 4.2.0)
# generics               0.1.2    2022-01-31 [2] CRAN (R 4.2.0)
# GenomeInfoDb         * 1.31.7   2022-04-01 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-02 [2] Bioconductor
# GenomicRanges        * 1.47.6   2022-01-12 [2] Bioconductor
# ggplot2              * 3.3.5    2021-06-25 [2] CRAN (R 4.2.0)
# glue                   1.6.2    2022-02-24 [2] CRAN (R 4.2.0)
# googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.2.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# haven                  2.4.3    2021-08-04 [2] CRAN (R 4.2.0)
# here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.2.0)
# hms                    1.1.1    2021-09-26 [2] CRAN (R 4.2.0)
# httr                   1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
# IRanges              * 2.29.1   2021-11-16 [2] Bioconductor
# jaffelab             * 0.99.32  2022-04-06 [1] Github (LieberInstitute/jaffelab@7b7afe3)
# jsonlite               1.8.0    2022-02-22 [2] CRAN (R 4.2.0)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.2.0)
# lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.2.0)
# limma                  3.51.5   2022-02-17 [2] Bioconductor
# lubridate              1.8.0    2021-10-07 [2] CRAN (R 4.2.0)
# magrittr               2.0.3    2022-03-30 [2] CRAN (R 4.2.0)
# MASS                   7.3-56   2022-03-23 [3] CRAN (R 4.2.0)
# Matrix                 1.4-1    2022-03-23 [3] CRAN (R 4.2.0)
# MatrixGenerics       * 1.7.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.2.0)
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.1.0)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# pillar                 1.7.0    2022-02-01 [2] CRAN (R 4.2.0)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R6                     2.5.1    2021-08-19 [2] CRAN (R 4.2.0)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.2.0)
# RColorBrewer           1.1-3    2022-04-03 [2] CRAN (R 4.2.0)
# RCurl                  1.98-1.6 2022-02-08 [2] CRAN (R 4.2.0)
# readr                * 2.1.2    2022-01-30 [2] CRAN (R 4.2.0)
# readxl                 1.4.0    2022-03-28 [2] CRAN (R 4.2.0)
# reprex                 2.0.1    2021-08-05 [2] CRAN (R 4.2.0)
# rlang                  1.0.2    2022-03-04 [2] CRAN (R 4.2.0)
# rprojroot              2.0.3    2022-04-02 [2] CRAN (R 4.2.0)
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.1.0)
# rvest                  1.0.2    2021-10-16 [2] CRAN (R 4.2.0)
# S4Vectors            * 0.33.16  2022-04-01 [2] Bioconductor
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# segmented              1.4-1    2022-03-24 [1] CRAN (R 4.2.0)
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.2.0)
# SingleCellExperiment * 1.17.2   2021-11-18 [2] Bioconductor
# stringi                1.7.6    2021-11-29 [2] CRAN (R 4.2.0)
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.25.3   2021-12-08 [2] Bioconductor
# tibble               * 3.1.6    2021-11-07 [2] CRAN (R 4.2.0)
# tidyr                * 1.2.0    2022-02-01 [2] CRAN (R 4.2.0)
# tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.2.0)
# tidyverse            * 1.3.1    2021-04-15 [2] CRAN (R 4.2.0)
# TREG                 * 0.99.4   2022-04-06 [1] Github (LieberInstitute/TREG@8f1dfd3)
# tzdb                   0.3.0    2022-03-28 [2] CRAN (R 4.2.0)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.2.0)
# vctrs                  0.4.0    2022-03-30 [2] CRAN (R 4.2.0)
# withr                  2.5.0    2022-03-03 [2] CRAN (R 4.2.0)
# xml2                   1.3.3    2021-11-30 [2] CRAN (R 4.2.0)
# XVector                0.35.0   2021-10-26 [2] Bioconductor
# zlibbioc               1.41.0   2021-10-26 [2] Bioconductor
#
# [1] /users/lhuuki/R/devel
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/R/devel/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/R/devel/lib64/R/library
#
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
