library("tidyverse")
library("SingleCellExperiment")
library("scater")
library("jaffelab")
library("VennDiagram")
library("ggrepel")
library("here")
library("sessioninfo")

plot_dir <- "plots/01_find_tregs"

#### Fig 1. Demo Rank Distribution ####
set.seed(10)
demo_rank_data <- tibble(
    sample = 1:100,
    `high\ninvariance` = round(rnorm(100, 75, 3), digits = 0),
    `low\ninvariance` = sample(1:100, 100, replace = TRUE)
) %>%
    pivot_longer(!sample, names_to = "Gene", values_to = "Expression Rank")

rank_demo_violin <- ggplot(demo_rank_data, aes(x = Gene, y = `Expression Rank`)) +
    geom_violin(fill = "light grey") +
    theme_bw() +
    theme(text = element_text(size = 15))

# ggsave(rank_demo_violin, filename = here(plot_dir, "main_pdf","fig1_rank_violin_demo.png"), width = 3, height = 6)
ggsave(rank_demo_violin, filename = here(plot_dir, "main_pdf", "fig1_rank_violin_demo.pdf"), width = 3, height = 3)

#### Prep sce data ####
load(here("raw-data", "sce_pan.v2.Rdata"), verbose = TRUE)
load(here("processed-data", "01_find_tregs", "rank_invar.Rdata"), verbose = TRUE)
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)

table(sce_pan$cellType.Broad)

## Capitalize Regions
sce_pan$region2 <- toupper(sce_pan$region)
sce_pan$region2[sce_pan$region2 == "SACC"] <- "sACC"
sce_pan$region2[sce_pan$region2 == "NAC"] <- "NAc"
## Rare cell types
sce_pan$region2[sce_pan$cellType.Broad %in% c("Macro", "Mural", "Endo", "Tcell")] <- "combined"
table(sce_pan$region2)

sce_pan$ctXregion <- paste0(sce_pan$cellType.Broad, "_", sce_pan$region2)

pd <- as.data.frame(colData(sce_pan))

#### Annotate Genes ####
gene_symbol <- rowData(sce_pan) %>%
    as.data.frame() %>%
    select(gene = gene_id, Symbol)

### TREGs Chosen from top 10 RI gene for probe availability
treg_list <- c("MALAT1", "AKT3", "ARID1B")

# classic HK genes
dotdotdot_genes <- c(
    "POLR2A", "PPIB", "UBC", "HPRT1", "ACTB", "TUBB3", "BIN1", "LDHA",
    "GAPDH", "PGK1", "BHLHE22", "CPLX2"
)

data_driven_HK <- c("NDUFB4", "NDUFB1", "GSTO1", "AMZ2", "POLR2I", "NDUFA3", "RRAGA", "POMP")

(candidate_genes <- gene_symbol %>%
    filter(Symbol %in% c(treg_list, "POLR2A")))


genes_of_interest <- tibble(
    Symbol = c(treg_list, dotdotdot_genes, data_driven_HK),
    gene_anno = c(
        rep("TREG Candidate", length(treg_list)),
        rep("Classic HK", length(dotdotdot_genes)),
        rep("Data Driven HK", length(data_driven_HK))
    )
) %>%
    left_join(gene_symbol)

write.csv(genes_of_interest, here("processed-data", "01_find_tregs", "supp_tables", "genes_of_interest.csv"))

#### Proportion zero plots ####
propZero_limit <- 0.75
## overall density

gene_propZero_long <- gene_propZero %>%
    rownames_to_column("gene") %>%
    pivot_longer(!gene, values_to = "propZero") %>%
    separate(name, into = c("cellType", "region"))

gene_propZero_long$cellType <- factor(gene_propZero_long$cellType, levels = levels(sce_pan$cellType.Broad))
gene_propZero_long$region <- factor(gene_propZero_long$region, levels = c("AMY", "DLPFC", "HPC", "NAc", "sACC", "combined"))

#### Fig 2 ####
propZero_density <- gene_propZero_long %>%
    filter(region != "combined") %>%
    ggplot(aes(x = propZero, fill = cellType)) +
    geom_histogram(binwidth = 0.05, color = "black", size = .2) +
    scale_fill_manual(values = cell_colors) +
    facet_grid(cellType ~ region) +
    labs(x = "Proportion Zero per Group", y = "Number of Genes") +
    geom_vline(xintercept = propZero_limit, color = "red", linetype = "dashed") +
    theme_bw() +
    theme(
        legend.position = "none",
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1)
    ) +
    scale_y_continuous(breaks = seq(0, 3000, 1500)) +
    scale_x_continuous(breaks = seq(0, 1, .5))

# ggsave(propZero_density, filename = here(plot_dir, "main_pdf","fig2_propZero_density.png"), width = 6, height = 7)
ggsave(propZero_density, filename = here(plot_dir, "main_pdf", "fig2_propZero_density.pdf"), width = 6, height = 7)


propZero_density_rare <- gene_propZero_long %>%
    filter(region == "combined") %>%
    ggplot(aes(x = propZero, fill = cellType)) +
    geom_histogram(binwidth = 0.05, color = "black", size = .2) +
    scale_fill_manual(values = cell_colors) +
    facet_wrap(~cellType, nrow = 1) +
    labs(x = "Proportion Zero per Group", y = "Number of Genes") +
    geom_vline(xintercept = propZero_limit, color = "red", linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(propZero_density_rare, filename = here(plot_dir, "supp_pdf", "fig_supp_propZero_density_rare_ct.png"), width = 7, height = 2)

## filter
# sce_pan <- sce_pan[genes_filtered,]
# (n_genes <- c(n_genes, "zero_filter" = nrow(sce_pan)))
# sce_pan$sum_counts_zero_filter <- colSums(assays(sce_pan)$counts)
# dim(sce_pan)
# # [1]   877 70527

# propZero_density_filtered <- gene_propZero_long %>%
#   filter(gene %in% genes_filtered) %>%
#   ggplot(aes(x = propZero, fill = cellType)) +
#   geom_histogram(binwidth = 0.05, color = "black",size=.2) +
#   scale_fill_manual(values = cell_colors) +
#   facet_grid(region~cellType) +
#   geom_vline(xintercept = propZero_limit, color = "red", linetype = "dashed") +
#   theme_bw()+
#   theme(legend.position = "none", te) +
#   labs( title = paste0("Proportion Zero Distribution - filter for ", propZero_limit, " max non-zero"),
#         x = "Proportion Zero per Group", y = "Number of Genes")
#
# ggsave(propZero_density_filtered, filename = here(plot_dir, "propZero_density_filter.png"), width = 10)

#### Demo filtering ####
filter_demo <- gene_symbol %>%
    filter(Symbol %in% c(treg_list, "POLR2A")) %>%
    left_join(gene_propZero_long)

filter_anno <- filter_demo %>%
    group_by(Symbol) %>%
    summarise(max = max(propZero))

filter_demo_scatter <- ggplot(filter_demo, aes(x = Symbol, y = propZero)) +
    geom_jitter(aes(color = cellType), width = 0.1) +
    scale_color_manual(values = cell_colors, name = "Cell Type") +
    labs(x = "Gene", y = "Proportion Zero") +
    geom_hline(yintercept = propZero_limit, color = "red", linetype = "dashed") +
    theme_bw() +
    theme(
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, face = "italic")
    )

# ggsave(filter_demo_scatter, filename = here(plot_dir, "main_pdf","fig2b_propZero_filter_demo.png"), width = 4, height = 7)
ggsave(filter_demo_scatter, filename = here(plot_dir, "main_pdf", "fig2b_propZero_filter_demo.pdf"), width = 4, height = 7)

#### Fig 3. Real Rank Distribution ####
load(here("processed-data", "01_find_tregs", "rank_df_subset.Rdata"), verbose = TRUE)
# 4 canidiate genes * 70k nuclei
dim(rank_df_subset)
# [1]     4 70527
corner(rank_df_subset)

rank_long <- rank_df_subset %>%
    rownames_to_column("gene") %>%
    pivot_longer(!gene, names_to = "Sample", values_to = "rank") %>%
    left_join(gene_symbol) %>%
    left_join(pd %>% select(Sample = uniqueID, cellType.Broad))

rank_violin <- ggplot(rank_long, aes(x = Symbol, y = rank)) +
    geom_violin(fill = "light grey", scale = "width") +
    labs(x = "Gene") +
    labs(y = "Expression Rank") +
    theme_bw() +
    theme(
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, face = "italic")
    )

# ggsave(rank_violin, filename = here(plot_dir, "main_pdf","fig3_rank_violin.png"), width = 3, height = 6)
ggsave(rank_violin, filename = here(plot_dir, "main_pdf", "fig3_rank_violin.pdf"), width = 3, height = 6)

rank_violin_ct <- ggplot(rank_long, aes(x = cellType.Broad, y = rank, fill = cellType.Broad)) +
    geom_violin(scale = "width") +
    labs(x = "Cell Type") +
    facet_wrap(~Symbol, ncol = 2) +
    scale_fill_manual(values = cell_colors) +
    labs(y = "Expression Rank") +
    theme_bw() +
    theme(
        text = element_text(size = 15),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(face = "italic")
    )

# ggsave(rank_violin_ct, filename = here(plot_dir, "main_pdf","fig3_rank_violin_ct.png"), width = 5.5, height = 6)
ggsave(rank_violin_ct, filename = here(plot_dir, "main_pdf", "fig3_rank_violin_ct.pdf"), width = 5.5, height = 6)


#### Fig 3. Do Expression trends over cell type track in our HK genes? ####
hk_counts <- log2(as.matrix(assays(sce_pan[genes_of_interest$gene, ])$counts) + 1)

counts_long <- as.data.frame(hk_counts) %>%
    rownames_to_column("gene") %>%
    pivot_longer(!gene, names_to = "ID", values_to = "logcount") %>%
    left_join(pd %>% rownames_to_column("ID") %>% select(ID, cellType.Broad, sum)) %>%
    mutate(logsum = log2(sum)) %>%
    left_join(genes_of_interest) %>%
    left_join(invar_t %>% select(gene = ensembl_id, t, adj.P.Val))

counts_long$Symbol <- factor(counts_long$Symbol)
counts_long$Symbol <- fct_reorder(counts_long$Symbol, counts_long$t, max)
levels(counts_long$Symbol)

model_anno <- invar_t %>%
    filter(gene_anno != "other") %>%
    mutate(anno = paste("t=", round(t, 1), "\nFDR=", scales::scientific(adj.P.Val, didgits = 3)))

model_anno$Symbol <- factor(model_anno$Symbol, levels = levels(counts_long$Symbol))

hk_sum_scatter <- counts_long %>%
    ggplot(aes(logsum, logcount)) +
    geom_point(aes(color = cellType.Broad), alpha = 0.5, size = 0.5) +
    facet_wrap(~Symbol) +
    geom_text(data = model_anno, aes(label = anno), x = 7.5, y = 12, vjust = "inward", hjust = "inward") +
    scale_color_manual(values = cell_colors) +
    theme_bw() +
    labs(x = "log2(sum)", y = "log2(count + 1)")

ggsave(hk_sum_scatter, filename = here(plot_dir, "explore", "hk_sum_scatter.png"), width = 12, height = 10)

hk_sum_smooth <- counts_long %>%
    ggplot() +
    geom_smooth(aes(logsum, logcount, color = cellType.Broad), method = "lm") +
    facet_wrap(~Symbol) +
    geom_text(data = model_anno, aes(label = anno), x = 7.5, y = 12, vjust = "inward", hjust = "inward") +
    scale_color_manual(values = cell_colors) +
    theme_bw() +
    labs(x = "log2(sum)", y = "log2(count + 1)")

ggsave(hk_sum_smooth, filename = here(plot_dir, "explore", "hk_sum_smooth.png"), width = 12, height = 10)

## Figure 3c
hk_sum_scatter_main <- counts_long %>%
    filter(Symbol == "MALAT1") %>%
    ggplot(aes(logsum, logcount, color = cellType.Broad)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_smooth(method = "lm") +
    scale_color_manual(values = cell_colors) +
    facet_wrap(~Symbol) +
    theme_bw() +
    coord_equal() +
    labs(x = "log2(sum)", y = "log2(count + 1)") +
    theme(
        legend.position = "None",
        text = element_text(size = 15),
        strip.text.x = element_text(face = "italic")
    )

# ggsave(hk_sum_scatter_main, filename = here("plots","01_find_tregs","main_pdf","fig3_hk_sum_scatter_MALAT1.png"), width = 4.5, height = 5)
ggsave(hk_sum_scatter_main, filename = here(plot_dir, "main_pdf", "fig3_hk_sum_scatter_MALAT1.pdf"), width = 4.5, height = 5)


hk_sum_smooth2 <- counts_long %>%
    filter(Symbol %in% c("AKT3", "ARID1B", "POLR2A")) %>%
    ggplot() +
    geom_smooth(aes(logsum, logcount, color = cellType.Broad), method = "lm") +
    facet_wrap(~Symbol, scales = "free_y", nrow = 1) +
    scale_color_manual(values = cell_colors) +
    theme_bw() +
    labs(x = "log2(sum)", y = "log2(count + 1)") +
    theme(legend.position = "None")

ggsave(hk_sum_smooth2, filename = here("plots", "01_find_tregs", "explore", "hk_sum_smooth2.png"), width = 12, height = 5)

## Figure 3d
hk_sum_smooth2_main <- counts_long %>%
    filter(Symbol %in% c("AKT3", "ARID1B", "POLR2A")) %>%
    ggplot() +
    geom_smooth(aes(logsum, logcount, color = cellType.Broad), method = "lm") +
    facet_wrap(~Symbol, scales = "free_y", ncol = 1) +
    scale_color_manual(values = cell_colors) +
    theme_bw() +
    labs(x = "log2(sum)", y = "log2(count + 1)") +
    theme(
        legend.position = "None", text = element_text(size = 15),
        strip.text.x = element_text(face = "italic")
    )

# ggsave(hk_sum_smooth2_main, filename = here("plots","01_find_tregs","main_pdf","fig3_hk_sum_smooth2.png"), width = 4, height = 5)
ggsave(hk_sum_smooth2_main, filename = here("plots", "01_find_tregs", "main_pdf", "fig3_hk_sum_smooth2.pdf"), width = 4, height = 5)

gene_slopes <- counts_long %>%
    group_by(Symbol) %>%
    do(fitQSV = broom::tidy(lm(logcount ~ logsum, data = .))) %>%
    unnest(fitQSV) %>%
    select(Symbol, term, estimate) %>%
    mutate(term = gsub("[^[:alnum:] ]", "", term)) %>%
    pivot_wider(Symbol, names_from = "term", values_from = "estimate") %>%
    column_to_rownames("Symbol")

hk_sum_scatter_AKT3 <- counts_long %>%
    filter(Symbol == "AKT3") %>%
    ggplot(aes(logsum, logcount, color = cellType.Broad)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_smooth(method = "lm", color = "black") +
    geom_abline(slope = 0.8672058, intercept = -8.7700637, color = "dark red") +
    facet_wrap(~cellType.Broad, scales = "free") +
    scale_color_manual(values = cell_colors) +
    theme_bw()

ggsave(hk_sum_scatter_AKT3, filename = here("plots", "01_find_tregs", "explore", "hk_sum_scatter_AKT3.png"), width = 11)

pdf(here("plots", "01_find_tregs", "explore", "hk_sum_scatter_smooth.pdf"), width = 11, height = 8)
for (g in levels(counts_long$Symbol)) {
    message(g)
    st <- model_anno %>%
        filter(Symbol == g) %>%
        pull(anno)

    m <- gene_slopes[g, "logsum"]
    b <- gene_slopes[g, "Intercept"]

    c <- counts_long %>%
        filter(Symbol == g) %>%
        ggplot(aes(logsum, logcount, color = cellType.Broad)) +
        geom_point(alpha = 0.5, size = 0.5) +
        geom_smooth(method = "lm", color = "black") +
        facet_wrap(~cellType.Broad, scales = "free") +
        geom_abline(slope = m, intercept = b, color = "dark red") +
        scale_color_manual(values = cell_colors) +
        theme_bw() +
        labs(title = g, subtitle = st, x = "log2(sum)", y = "log2(count + 1)")

    print(c)
}
dev.off()


# sgejobs::job_single('TREG_snRNA-seq_plots', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript TREG_snRNA-seq_plots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

# [1] "Reproducibility information:"
# [1] "2022-03-03 17:04:57 EST"
# user  system elapsed
# 304.203  15.204 325.654
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
# date     2022-03-03
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# backports              1.4.1    2021-12-13 [2] CRAN (R 4.1.2)
# beachmat               2.10.0   2021-10-26 [2] Bioconductor
# beeswarm               0.4.0    2021-06-01 [2] CRAN (R 4.1.2)
# Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
# BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
# BiocNeighbors          1.12.0   2021-10-26 [2] Bioconductor
# BiocParallel           1.28.3   2021-12-09 [2] Bioconductor
# BiocSingular           1.10.0   2021-10-26 [2] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# broom                  0.7.12   2022-01-28 [2] CRAN (R 4.1.2)
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.1.0)
# cli                    3.2.0    2022-02-14 [2] CRAN (R 4.1.2)
# colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.1.2)
# crayon                 1.5.0    2022-02-14 [2] CRAN (R 4.1.2)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
# dbplyr                 2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
# digest                 0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
# dplyr                * 1.0.8    2022-02-08 [2] CRAN (R 4.1.2)
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  1.0.2    2022-01-14 [2] CRAN (R 4.1.2)
# farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
# forcats              * 0.5.1    2021-01-27 [2] CRAN (R 4.1.0)
# formatR                1.11     2021-06-01 [2] CRAN (R 4.1.0)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
# futile.logger        * 1.4.3    2016-07-10 [2] CRAN (R 4.1.0)
# futile.options         1.0.1    2018-04-20 [2] CRAN (R 4.1.0)
# gargle                 1.2.0    2021-07-02 [2] CRAN (R 4.1.0)
# generics               0.1.2    2022-01-31 [2] CRAN (R 4.1.2)
# GenomeInfoDb         * 1.30.1   2022-01-30 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
# GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
# ggbeeswarm             0.6.0    2017-08-07 [2] CRAN (R 4.1.2)
# ggplot2              * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
# ggrepel              * 0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
# glue                   1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
# googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
# gridExtra              2.3      2017-09-09 [2] CRAN (R 4.1.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# haven                  2.4.3    2021-08-04 [2] CRAN (R 4.1.1)
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# hms                    1.1.1    2021-09-26 [2] CRAN (R 4.1.2)
# httr                   1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# irlba                  2.3.5    2021-12-06 [1] CRAN (R 4.1.2)
# jaffelab             * 0.99.31  2021-10-29 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# jsonlite               1.8.0    2022-02-22 [2] CRAN (R 4.1.2)
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
# lambda.r               1.2.4    2019-09-18 [2] CRAN (R 4.1.0)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [1] CRAN (R 4.1.1)
# limma                  3.50.1   2022-02-17 [2] Bioconductor
# lubridate              1.8.0    2021-10-07 [2] CRAN (R 4.1.2)
# magrittr               2.0.2    2022-01-26 [2] CRAN (R 4.1.2)
# MASS                   7.3-55   2022-01-13 [3] CRAN (R 4.1.2)
# Matrix                 1.4-0    2021-12-08 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
# mgcv                   1.8-39   2022-02-24 [3] CRAN (R 4.1.2)
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.1.0)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# nlme                   3.1-155  2022-01-13 [3] CRAN (R 4.1.2)
# pillar                 1.7.0    2022-02-01 [1] CRAN (R 4.1.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R6                     2.5.1    2021-08-19 [1] CRAN (R 4.1.1)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.1)
# ragg                   1.2.2    2022-02-21 [2] CRAN (R 4.1.2)
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.8    2022-01-13 [2] CRAN (R 4.1.2)
# RCurl                  1.98-1.6 2022-02-08 [2] CRAN (R 4.1.2)
# readr                * 2.1.2    2022-01-30 [2] CRAN (R 4.1.2)
# readxl                 1.3.1    2019-03-13 [2] CRAN (R 4.1.0)
# reprex                 2.0.1    2021-08-05 [2] CRAN (R 4.1.1)
# rlang                  1.0.1    2022-02-03 [2] CRAN (R 4.1.2)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.1.0)
# rsvd                   1.0.5    2021-04-16 [2] CRAN (R 4.1.2)
# rvest                  1.0.2    2021-10-16 [2] CRAN (R 4.1.2)
# S4Vectors            * 0.32.3   2021-11-21 [2] Bioconductor
# ScaledMatrix           1.2.0    2021-10-26 [2] Bioconductor
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# scater               * 1.22.0   2021-10-26 [2] Bioconductor
# scuttle              * 1.4.0    2021-10-26 [2] Bioconductor
# segmented              1.4-0    2022-01-28 [1] CRAN (R 4.1.2)
# sessioninfo            1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
# SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
# sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
# stringi                1.7.6    2021-11-29 [2] CRAN (R 4.1.2)
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
# systemfonts            1.0.4    2022-02-11 [2] CRAN (R 4.1.2)
# textshaping            0.3.6    2021-10-13 [2] CRAN (R 4.1.2)
# tibble               * 3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
# tidyr                * 1.2.0    2022-02-01 [2] CRAN (R 4.1.2)
# tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
# tidyverse            * 1.3.1    2021-04-15 [2] CRAN (R 4.1.0)
# tzdb                   0.2.0    2021-10-27 [2] CRAN (R 4.1.2)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# VennDiagram          * 1.7.1    2021-12-02 [2] CRAN (R 4.1.2)
# vipor                  0.4.5    2017-03-22 [2] CRAN (R 4.1.2)
# viridis                0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
# viridisLite            0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
# withr                  2.4.3    2021-11-30 [2] CRAN (R 4.1.2)
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
#   Thu Mar  3 17:05:01 EST 2022
