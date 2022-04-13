library("tidyverse")
library("here")
library("sessioninfo")
library("ggrepel")
library("UpSetR")
library("DeconvoBuddies")
library("jaffelab")
library("org.Hs.eg.db")
library("clusterProfiler")

genes_of_interest <- read.csv(here("processed-data", "01_find_tregs", "supp_tables", "genes_of_interest.csv"), row.names = 1)
gene_metrics <- read.csv(here("processed-data", "01_find_tregs", "supp_tables", "gene_metrics.csv"), row.names = 1)
head(gene_metrics)

## annotate gene metrics
gene_metrics <- gene_metrics %>%
    left_join(genes_of_interest %>% select(Symbol, `Gene Type` = gene_anno))

gene_metrics %>%
    filter(!is.na(`Gene Type`))


## Eisenberg et al. HKGs
eb_HKG <- c("C1orf43", "CHMP2A", "EMC7", "GPI", "PSMB2", "PSMB4", "RAB7A", "REEP5", "SNRPD3", "VCP", "VPS29")
# Symbol      ensembl_id top50 max_PropZero PropZero_filter rank_invar Gene Type
#     1    PSMB2 ENSG00000126067  TRUE    1.0000000           FALSE         NA      <NA>
#     2    PSMB4 ENSG00000159377  TRUE    0.9793814           FALSE         NA      <NA>
#     3  C1orf43 ENSG00000143612  TRUE    0.9062500           FALSE         NA      <NA>
#     4    RAB7A ENSG00000075785  TRUE    0.6451613            TRUE        695      <NA>
#     5    REEP5 ENSG00000129625  TRUE    0.8750000           FALSE         NA      <NA>
#     6      VCP ENSG00000165280  TRUE    0.9375000           FALSE         NA      <NA>
#     7    VPS29 ENSG00000111237  TRUE    0.8636364           FALSE         NA      <NA>
#     8     EMC7 ENSG00000134153  TRUE    0.9677419           FALSE         NA      <NA>
#     9      GPI ENSG00000105220  TRUE    0.8939394           FALSE         NA      <NA>
#     10  CHMP2A ENSG00000130724 FALSE           NA           FALSE         NA      <NA>
#     11  SNRPD3 ENSG00000100028  TRUE    0.9298469           FALSE         NA      <NA>


gene_metrics %>% filter(Symbol %in% eb_HKG)

## Add T-stats
load(here("processed-data", "01_find_tregs", "lmfit.Rdata"), verbose = TRUE) # tt
tt2 <- tt %>%
    select(t) %>%
    rownames_to_column("ensembl_id")

gene_metrics2 <- gene_metrics %>%
    left_join(tt2)

## export version
## Exclude marker annotation for now
head(gene_metrics2)
#       Symbol      ensembl_id top50 max_PropZero PropZero_filter rank_invar Gene Type          t
# 1 AL627309.1 ENSG00000238009  TRUE      1.00000           FALSE         NA      <NA>  51.188405
# 2 AC114498.1 ENSG00000235146 FALSE           NA           FALSE         NA      <NA>   9.844346
# 3 AL669831.5 ENSG00000237491  TRUE      0.96875           FALSE         NA      <NA> 106.687296
# 4  LINC00115 ENSG00000225880 FALSE           NA           FALSE         NA      <NA>  18.115699
# 5     FAM41C ENSG00000230368 FALSE           NA           FALSE         NA      <NA>   8.931643
# 6 AL645608.7 ENSG00000272438 FALSE           NA           FALSE         NA      <NA>   6.048259
write.csv(gene_metrics2, here("processed-data", "01_find_tregs", "supp_tables", "gene_metrics2.csv"))


#### t-stat plots ####
gene_metrics3 <- gene_metrics2 %>%
    replace_na(list(`Gene Type` = "None")) %>%
    mutate(
        gene_anno = case_when(
            !top50 ~ "Fail 50% Exp.",
            !PropZero_filter ~ "Fail Prop. Zero",
            TRUE ~ "Evaluated RI"
        ),
        label1 = Symbol %in% c("ARID1B", "AKT3", "MALAT1", "POLR2A"),
        label2 = `Gene Type` != "None",
        alpha = `Gene Type` != "None",
        `Gene Type` = factor(`Gene Type`, levels = c("TREG Candidate", "Classic HK", "Data Driven HK", "None"))
    )

gene_metrics3 %>% dplyr::count(gene_anno, `Gene Type`)
#         gene_anno      Gene Type     n
# 1    Evaluated RI  TREG Candidate     3
# 2    Evaluated RI           None   874
# 3   Fail 50% Exp.     Classic HK     2
# 4   Fail 50% Exp. Data Driven HK     1
# 5   Fail 50% Exp.           None 11516
# 6 Fail Prop. Zero     Classic HK    10
# 7 Fail Prop. Zero Data Driven HK     7
# 8 Fail Prop. Zero           None 10625

## t-stat plots
plot_dir <- "plots/01_find_tregs"
gene_type_colors <- create_cell_colors(levels(gene_metrics3$`Gene Type`), pallet = "gg")
gene_type_colors["None"] <- "grey"
gene_type_colors["TRUE"] <- "blue"

invar_t_scatter <- gene_metrics3 %>%
    filter(gene_anno == "Evaluated RI") %>%
    ggplot(aes(x = rank_invar, y = t)) +
    geom_point(alpha = 0.2, color = "grey") +
    geom_point(
        data = filter(gene_metrics3, label2, gene_anno == "Evaluated RI"),
        aes(color = `Gene Type`)
    ) +
    geom_point(data = filter(gene_metrics3, label1, gene_anno == "Evaluated RI"), shape = 21, color = "black") +
    geom_text_repel(aes(
        label = ifelse(label2, paste0("italic('", Symbol, "')"), NA),
        color = label1
    ),
    size = 3, parse = TRUE, show.legend = FALSE
    ) +
    scale_color_manual(values = gene_type_colors, drop = FALSE) +
    ylim(0, 400) +
    theme_bw() +
    facet_wrap(~gene_anno) +
    labs(x = "Rank Invariance", y = "Total RNA t-statistic") +
    theme(text = element_text(size = 15))

ggsave(invar_t_scatter, filename = here(plot_dir, "explore", "rank_invar_t_scatter.png"), width = 6)
ggsave(invar_t_scatter, filename = here(plot_dir, "supp_pdf", "rank_invar_t_scatter.pdf"), width = 6)

## Denisity/jitter plot
pos <- position_jitter(width = 0.3, seed = 2)
invar_t_density <- gene_metrics3 %>%
    filter(gene_anno != "Evaluated RI") %>%
    ggplot(aes(x = `Gene Type`, y = t)) +
    geom_point(aes(color = `Gene Type`, alpha = alpha), position = pos) +
    # geom_point(aes(alpha = alpha, color = `Gene Type`), position = pos) +
    # geom_point(data = filter(gene_metrics3, label2, gene_anno != "Evaluated RI"),
    #            aes(fill= `Gene Type`), shape=21, color = "black", position = pos)+
    geom_text_repel(aes(
        label = ifelse(label2, paste0("italic('", Symbol, "')"), NA),
        color = label1
    ),
    size = 3, position = pos, parse = TRUE
    ) +
    scale_color_manual(values = gene_type_colors, drop = FALSE) +
    facet_wrap(~gene_anno) +
    ylim(0, 400) +
    theme_bw() +
    labs(y = "Total RNA t-statistic") +
    theme(
        legend.position = "none",
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()
    )

ggsave(invar_t_density, filename = here(plot_dir, "explore", "rank_invar_t_densitiy.png"), width = 6, height = 7.5)
ggsave(invar_t_density, filename = here(plot_dir, "supp_pdf", "rank_invar_t_densitiy.pdf"), width = 6, height = 7.5)

#### GO Enrichment ####

## define universe
all_entrez <- bitr(gene_metrics2$ensembl_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
nrow(all_entrez)
# [1] 18296
u <- all_entrez$ENTREZ

## List top 10% (87) of the rank invar genes
top_RI <- gene_metrics2 %>%
    filter(!is.na(rank_invar)) %>%
    arrange(-rank_invar) %>%
    head(sum(!is.na(gene_metrics2$rank_invar)) %/% 10) %>%
    pull(ensembl_id)

top_RI_entrez <- bitr(top_RI, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZ
length(top_RI_entrez)
# [1] 86

## Run Enrichment
ont <- c("BP", "CC", "MF")
names(ont) <- ont
go_output <- map(ont, ~ compareCluster(
    geneClusters = list(topRI = top_RI_entrez),
    univ = u,
    OrgDb = "org.Hs.eg.db",
    fun = "enrichGO",
    ont = .x
))

pdf(here(plot_dir, "explore", "go_all_top87.pdf"))
map2(go_output, names(go_output), ~ dotplot(.x, title = .y))
dev.off()

save(go_output, file = here("processed-data", "01_find_tregs", "enrichGO.Rdata"))

go_df <- map2(go_output, names(go_output), ~ as.data.frame(.x) %>% mutate(ont = .y))
go_df_all <- do.call("rbind", go_df) %>%
    as_tibble() %>%
    separate(GeneRatio, into = c("n1", "n2"), sep = "/") %>%
    mutate(
        GeneRatio = Count / as.integer(n2),
        Cluster = paste0(Cluster, "\n(", n2, ")"),
        Description = gsub("\n$", "", gsub("(.{1,20})(\\s|$)", "\\1\n", Description))
    )

gg_dotplot <- go_df_all %>%
    ggplot(aes(x = Cluster, y = Description, color = p.adjust, size = GeneRatio)) +
    geom_point() +
    facet_wrap(~ont, scales = "free") +
    scale_colour_gradient(low = "red", high = "blue") +
    theme_bw()

ggsave(gg_dotplot, filename = here(plot_dir, "explore", "go_gg_top87.png"), width = 11)
ggsave(gg_dotplot, filename = here(plot_dir, "supp_pdf", "go_gg_top87.pdf"), width = 11)

#### Marker Test ####
## Marker genes from Tran Maynard et al.
## https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/810b47364af4c8afe426bd2a6b559bd6a9f1cc98/tables/revision/top40genesLists_DLPFC-n3_cellType_SN-LEVEL-tests_LAH2020.csv
tran_markers <- read.csv(here("raw-data", "top40genesLists_DLPFC-n3_cellType_SN-LEVEL-tests_LAH2020.csv"))

tran_markers_long <- tran_markers %>%
    rownames_to_column("rank") %>%
    pivot_longer(cols = !rank, names_to = "Marker", values_to = "Symbol") %>%
    filter(Symbol != "") %>%
    separate(Marker, into = c("cellType", "test"), "_", extra = "merge", remove = FALSE) %>%
    separate(test, into = c("cellType_suf", "test"), "_", fill = "left") %>%
    mutate(
        cellType = ifelse(is.na(cellType_suf), cellType, paste0(cellType, "_", cellType_suf)),
        test = paste0("Tran_", test)
    ) %>%
    select(-cellType_suf)

tran_markers_long %>% count(cellType, test)

## Make table for Park markers
park_markers_list <- readLines(here("raw-data", "Park_markers_n69.txt"))
park_markers <- tibble(
    Symbol = park_markers_list,
    test = "Park"
)

markers_anno <- tran_markers_long %>%
    select(Symbol, test) %>%
    rbind(park_markers) %>%
    group_by(Symbol) %>%
    mutate(test = ordered(test, levels = c("Tran_1vAll", "Tran_pw", "Park"))) %>%
    arrange(test) %>%
    summarise(Marker = paste0(unique(test), collapse = ","))

# Marker                      n
# <chr>                   <int>
# 1 Park                       41
# 2 Tran_1vAll                203
# 3 Tran_1vAll,Park             8
# 4 Tran_1vAll,Tran_pw        396
# 5 Tran_1vAll,Tran_pw,Park    20
# 6 Tran_pw                   188

gene_metrics2 <- gene_metrics2 %>%
    left_join(markers_anno)

## plotting version
# gene_metrics_marker <- gene_metrics2 %>%
#     rename(gene_type = `Gene Type`) %>%
#     mutate(`Gene Type` = factor(ifelse(!is.na(Marker), Marker, gene_type),
#         levels = c(
#             "TREG Candidate", "Classic HK", "Data Driven HK",
#             "Tran_1vAll", "Tran_pw", "Park", "Tran_1vAll,Tran_pw",
#             "Tran_1vAll,Park", "Tran_1vAll,Tran_pw,Park", "None"
#         )
#     )) %>%
#     replace_na(list(`Gene Type` = "None"))

## no overlapping Marker + gene_type annotations
gene_metrics2 %>%
    count(`Gene Type`) %>%
    arrange(`Gene Type`)

## 11 genes are markers + have RI vals
gene_metrics_marker %>%
    filter(!is.na(Marker) & !is.na(rank_invar)) %>%
    arrange(-rank_invar)

#     Symbol      ensembl_id top50 max_PropZero PropZero_filter rank_invar gene_type         t             Marker
# 1     SIK3 ENSG00000160584  TRUE    0.3762887            TRUE      838.0      <NA> 288.31400 Tran_1vAll,Tran_pw * need to check this one out!
# 2      QKI ENSG00000112531  TRUE    0.4242424            TRUE      817.0      <NA> 197.21176            Tran_pw
# 3  TMEM165 ENSG00000134851  TRUE    0.6662404            TRUE      765.0      <NA> 172.69393 Tran_1vAll,Tran_pw
# 4    PCDH9 ENSG00000184226  TRUE    0.6875000            TRUE      675.0      <NA> 219.19276            Tran_pw
# 5     DLG1 ENSG00000075711  TRUE    0.4393939            TRUE      490.0      <NA> 199.71697            Tran_pw
# 6     CBLB ENSG00000114423  TRUE    0.7216495            TRUE      440.0      <NA> 118.95549            Tran_pw
# 7  PIP4K2A ENSG00000150867  TRUE    0.5806452            TRUE      283.0      <NA> 180.39774         Tran_1vAll
# 8    RAP1B ENSG00000127314  TRUE    0.7183682            TRUE      249.0      <NA> 107.27441            Tran_pw
# 9    PDE8A ENSG00000073417  TRUE    0.6818182            TRUE      224.0      <NA> 179.72187         Tran_1vAll
# 10    FRYL ENSG00000075539  TRUE    0.4232737            TRUE      136.5      <NA> 206.07683         Tran_1vAll
# 11    TLE4 ENSG00000106829  TRUE    0.6939891            TRUE       83.0      <NA>  99.10509            Tran_pw



#### Upset Plots ####
gene_metrics_filter <- gene_metrics %>% filter(!is.na(`Gene Type`))

gene_lists <- c(
    list(Park = park_markers_list),
    map(splitit(tran_markers_long$test), ~ tran_markers_long$Symbol[.x]),
    map(splitit(gene_metrics_filter$`Gene Type`), ~ gene_metrics_filter$Symbol[.x])
)

map(gene_lists, head)

pdf(here(plot_dir, "supp_pdf", "upset.pdf"))
upset(fromList(gene_lists), order.by = "freq", nset = length(gene_lists))
dev.off()


## Mean ratio markers
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/marker_stats_pan.v2.Rdata", verbose = TRUE)

marker_stats <- marker_stats %>%
    filter(rank_ratio <= 25)


gene_metrics2 %>%
    left_join(marker_stats %>% select(Symbol, cellType)) %>%
    count(`Gene Type`, cellType, PropZero_filter, !is.na(rank_invar))

#         Gene Type cellType !is.na(rank_invar)     n
# 1      Classic HK     <NA>              FALSE    12
# 2  Data Driven HK     <NA>              FALSE     8
# 3   TREG Canidate     <NA>               TRUE     3
# 4            <NA>    Astro              FALSE    11
# 5            <NA>     Endo              FALSE    25
# 6            <NA>    Macro              FALSE    23
# 7            <NA>    Micro              FALSE    32
# 8            <NA>    Mural              FALSE    48
# 9            <NA>    Oligo              FALSE     8
# 10           <NA>      OPC              FALSE    11
# 11           <NA>    Tcell              FALSE    11
# 12           <NA>    Excit              FALSE    46
# 13           <NA>    Inhib              FALSE    35
# 14           <NA>     <NA>              FALSE 21891
# 15           <NA>     <NA>               TRUE   874

# sgejobs::job_single('gene_check', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript gene_check.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# * 0.40.0   2021-10-26 [2] Bioconductor
# BiocNeighbors          1.12.0   2021-10-26 [2] Bioconductor
# BiocParallel           1.28.3   2021-12-09 [2] Bioconductor
# BiocSingular           1.10.0   2021-10-26 [2] Bioconductor
# Biostrings             2.62.0   2021-10-26 [2] Bioconductor
# bit                    4.0.4    2020-08-04 [2] CRAN (R 4.1.0)
# bit64                  4.0.5    2020-08-30 [2] CRAN (R 4.1.0)
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# blob                   1.2.2    2021-07-23 [2] CRAN (R 4.1.0)
# bluster                1.4.0    2021-10-26 [2] Bioconductor
# broom                  0.7.12   2022-01-28 [2] CRAN (R 4.1.2)
# cachem                 1.0.6    2021-08-19 [2] CRAN (R 4.1.2)
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.1.0)
# cli                    3.2.0    2022-02-14 [2] CRAN (R 4.1.2)
# cluster                2.1.3    2022-03-28 [3] CRAN (R 4.1.2)
# clusterProfiler      * 4.2.2    2022-01-13 [1] Bioconductor
# colorout             * 1.2-2    2021-12-03 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.1.2)
# crayon                 1.5.1    2022-03-26 [2] CRAN (R 4.1.2)
# data.table             1.14.2   2021-09-27 [2] CRAN (R 4.1.2)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
# dbplyr                 2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
# DeconvoBuddies       * 0.99.0   2021-12-03 [1] Github (lahuuki/DeconvoBuddies@1f6f074)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
# digest                 0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
# DO.db                  2.9      2021-12-03 [1] Bioconductor
# DOSE                   3.20.1   2021-11-18 [1] Bioconductor
# downloader             0.4      2015-07-09 [2] CRAN (R 4.1.0)
# dplyr                * 1.0.8    2022-02-08 [2] CRAN (R 4.1.2)
# dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.1.2)
# edgeR                  3.36.0   2021-10-26 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# enrichplot             1.14.1   2021-10-31 [1] Bioconductor
# fansi                  1.0.3    2022-03-24 [2] CRAN (R 4.1.2)
# farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
# fastmap                1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
# fastmatch              1.1-3    2021-07-23 [1] CRAN (R 4.1.2)
# fgsea                  1.20.0   2021-10-26 [1] Bioconductor
# forcats              * 0.5.1    2021-01-27 [2] CRAN (R 4.1.0)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
# gargle                 1.2.0    2021-07-02 [2] CRAN (R 4.1.0)
# generics               0.1.2    2022-01-31 [2] CRAN (R 4.1.2)
# GenomeInfoDb           1.30.1   2022-01-30 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
# GenomicRanges          1.46.1   2021-11-18 [2] Bioconductor
# ggforce                0.3.3    2021-03-05 [2] CRAN (R 4.1.0)
# ggfun                  0.0.5    2022-01-20 [1] CRAN (R 4.1.2)
# ggplot2              * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
# ggplotify              0.1.0    2021-09-02 [1] CRAN (R 4.1.2)
# ggraph                 2.0.5    2021-02-23 [2] CRAN (R 4.1.0)
# ggrepel              * 0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
# ggtree                 3.2.1    2021-11-16 [1] Bioconductor
# glue                   1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
# GO.db                  3.14.0   2021-11-01 [2] Bioconductor
# googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
# GOSemSim               2.20.0   2021-10-26 [1] Bioconductor
# graphlayouts           0.8.0    2022-01-03 [2] CRAN (R 4.1.2)
# gridExtra              2.3      2017-09-09 [2] CRAN (R 4.1.0)
# gridGraphics           0.5-1    2020-12-13 [1] CRAN (R 4.1.2)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# haven                  2.4.3    2021-08-04 [2] CRAN (R 4.1.1)
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# hms                    1.1.1    2021-09-26 [2] CRAN (R 4.1.2)
# httr                   1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
# igraph                 1.2.11   2022-01-04 [2] CRAN (R 4.1.2)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# irlba                  2.3.5    2021-12-06 [1] CRAN (R 4.1.2)
# jaffelab             * 0.99.31  2021-10-29 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# jsonlite               1.8.0    2022-02-22 [2] CRAN (R 4.1.2)
# KEGGREST               1.34.0   2021-10-26 [2] Bioconductor
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lazyeval               0.2.2    2019-03-15 [2] CRAN (R 4.1.0)
# lifecycle              1.0.1    2021-09-24 [1] CRAN (R 4.1.1)
# limma                  3.50.1   2022-02-17 [2] Bioconductor
# locfit                 1.5-9.5  2022-03-03 [2] CRAN (R 4.1.2)
# lubridate              1.8.0    2021-10-07 [2] CRAN (R 4.1.2)
# magrittr               2.0.2    2022-01-26 [2] CRAN (R 4.1.2)
# MASS                   7.3-56   2022-03-23 [3] CRAN (R 4.1.2)
# Matrix                 1.4-1    2022-03-23 [3] CRAN (R 4.1.2)
# MatrixGenerics         1.6.0    2021-10-26 [2] Bioconductor
# matrixStats            0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
# memoise                2.0.1    2021-11-26 [2] CRAN (R 4.1.2)
# metapod                1.2.0    2021-10-26 [2] Bioconductor
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.1.0)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# nlme                   3.1-157  2022-03-25 [3] CRAN (R 4.1.2)
# org.Hs.eg.db         * 3.14.0   2021-11-01 [2] Bioconductor
# patchwork              1.1.1    2020-12-17 [2] CRAN (R 4.1.2)
# pillar                 1.7.0    2022-02-01 [1] CRAN (R 4.1.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# plyr                   1.8.7    2022-03-24 [2] CRAN (R 4.1.2)
# png                    0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
# polyclip               1.10-0   2019-03-14 [2] CRAN (R 4.1.0)
# purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# qvalue                 2.26.0   2021-10-26 [2] Bioconductor
# R6                     2.5.1    2021-08-19 [1] CRAN (R 4.1.1)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.1)
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.8.3  2022-03-17 [2] CRAN (R 4.1.2)
# RCurl                  1.98-1.6 2022-02-08 [2] CRAN (R 4.1.2)
# readr                * 2.1.2    2022-01-30 [2] CRAN (R 4.1.2)
# readxl                 1.4.0    2022-03-28 [2] CRAN (R 4.1.2)
# reprex                 2.0.1    2021-08-05 [2] CRAN (R 4.1.1)
# reshape2               1.4.4    2020-04-09 [2] CRAN (R 4.1.0)
# rlang                  1.0.2    2022-03-04 [2] CRAN (R 4.1.2)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# RSQLite                2.2.11   2022-03-23 [2] CRAN (R 4.1.2)
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.1.0)
# rsvd                   1.0.5    2021-04-16 [2] CRAN (R 4.1.2)
# rvest                  1.0.2    2021-10-16 [2] CRAN (R 4.1.2)
# S4Vectors            * 0.32.4   2022-03-24 [2] Bioconductor
# ScaledMatrix           1.2.0    2021-10-26 [2] Bioconductor
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# scatterpie             0.1.7    2021-08-20 [1] CRAN (R 4.1.2)
# scran                  1.22.1   2021-11-14 [2] Bioconductor
# scuttle                1.4.0    2021-10-26 [2] Bioconductor
# segmented              1.4-0    2022-01-28 [1] CRAN (R 4.1.2)
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
# shadowtext             0.1.1    2022-01-10 [1] CRAN (R 4.1.2)
# SingleCellExperiment   1.16.0   2021-10-26 [2] Bioconductor
# sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
# statmod                1.4.36   2021-05-10 [2] CRAN (R 4.1.0)
# stringi                1.7.6    2021-11-29 [2] CRAN (R 4.1.2)
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment   1.24.0   2021-10-26 [2] Bioconductor
# tibble               * 3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
# tidygraph              1.2.0    2020-05-12 [2] CRAN (R 4.1.0)
# tidyr                * 1.2.0    2022-02-01 [2] CRAN (R 4.1.2)
# tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
# tidytree               0.3.8    2022-02-17 [1] CRAN (R 4.1.2)
# tidyverse            * 1.3.1    2021-04-15 [2] CRAN (R 4.1.0)
# treeio                 1.18.1   2021-11-14 [1] Bioconductor
# tweenr                 1.0.2    2021-03-23 [2] CRAN (R 4.1.0)
# tzdb                   0.3.0    2022-03-28 [2] CRAN (R 4.1.2)
# UpSetR               * 1.4.0    2019-05-22 [2] CRAN (R 4.1.2)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# viridis                0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
# viridisLite            0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
# withr                  2.5.0    2022-03-03 [2] CRAN (R 4.1.2)
# xml2                   1.3.3    2021-11-30 [2] CRAN (R 4.1.2)
# XVector                0.34.0   2021-10-26 [2] Bioconductor
# yulab.utils            0.0.4    2021-10-09 [1] CRAN (R 4.1.2)
# zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
#
# [1] /users/lhuuki/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
