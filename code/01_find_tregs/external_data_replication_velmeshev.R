library("SingleCellExperiment")
library("GGally")
library("tidyverse")
library("here")
library("TREG")
## Load Data
gene_metrics <- read.csv(here("processed-data", "01_find_tregs", "supp_tables", "gene_metrics2.csv"), row.names = 1)
cell_type_match <- read.csv(here("processed-data", "01_find_tregs", "Velmeshev_broad.csv")) %>%
    select(-Note)

load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/SCE_asd-velmeshev-etal_controls_multiBatchNorm.Rdata", verbose = TRUE)
dim(sce.asd)
# 36501 52556

## Populate rowData
rowData(sce.asd)

rowData(sce.asd) <- DataFrame(ensmblID = rownames(sce.asd))

#### Filter to Prop Zero Genes from 10x data ####
gene_set <- gene_metrics$ensembl_id[gene_metrics$PropZero_filter]
length(gene_set)
# [1] 877
table(gene_set %in% row.names(sce.asd))
# FALSE  TRUE
#     1   876
sce.asd <- sce.asd[gene_set[gene_set %in% row.names(sce.asd)], ]
dim(sce.asd)


#### Match Up Cell Types ####
pd <- colData(sce.asd) %>%
    as.data.frame() %>%
    left_join(cell_type_match)
colData(sce.asd) <- DataFrame(pd)

cat(levels(sce.asd$cluster), "\n")
table(sce.asd$cluster, sce.asd$cellType.Broad)
# Astro Endo Excit Inhib Micro Oligo  OPC
# AST-FB            2114    0     0     0     0     0    0
# AST-PP            2530    0     0     0     0     0    0
# Endothelial          0 1548     0     0     0     0    0
# IN-PV                0    0     0  2062     0     0    0
# IN-SST               0    0     0  2449     0     0    0
# IN-SV2C              0    0     0   887     0     0    0
# IN-VIP               0    0     0  2686     0     0    0
# L2/3                 0    0  6085     0     0     0    0
# L4                   0    0  3258     0     0     0    0
# L5/6                 0    0  1780     0     0     0    0
# L5/6-CC              0    0  2210     0     0     0    0
# Microglia            0    0     0     0  1822     0    0
# Neu-mat              0    0  1981     0     0     0    0
# Neu-NRGN-I           0    0  2027     0     0     0    0
# Neu-NRGN-II          0    0  4731     0     0     0    0
# Oligodendrocytes     0    0     0     0     0  9603    0
# OPC                  0    0     0     0     0     0 4783

summary(sce.asd$age)

## Run RI with all data, Velm clusters
rank_invar_velm <- rank_invariance_express(sce.asd, group_col = "cluster")
gene_metrics_velm <- as.data.frame(rank_invar_velm) %>%
    rownames_to_column("ensembl_id")

## Run RI with all data, broad cell types
rank_invar_velm_broad <- rank_invariance_express(sce.asd, group_col = "cellType.Broad")
gene_metrics_velm_broad <- as.data.frame(rank_invar_velm_broad) %>%
    rownames_to_column("ensembl_id")

## Run RI with 17+ data, broad cell types
sce.asd <- sce.asd[, sce.asd$age > 17]
dim(sce.asd)
# [1]   876 17784

rank_invar_velm_age <- rank_invariance_express(sce.asd, group_col = "cluster")
gene_metrics_velm_age <- as.data.frame(rank_invar_velm_age) %>%
    rownames_to_column("ensembl_id")

## Add to gene  metrics
gene_metrics <- gene_metrics %>%
    left_join(gene_metrics_velm) %>%
    left_join(gene_metrics_velm_broad) %>%
    left_join(gene_metrics_velm_age)

write.csv(gene_metrics, file = here("processed-data", "01_find_tregs", "gene_metrics_velm.csv"))

gene_metrics %>%
    filter(rank_invar > 867) %>%
    arrange(-rank_invar) %>%
    select(Symbol, starts_with("rank_invar"))

gene_metrics %>%
    filter(rank_invar < 10) %>%
    arrange(-rank_invar) %>%
    select(Symbol, starts_with("rank_invar"))

#### Plots ####
plot_dir <- "plots/01_find_tregs"

# scatter_RI <- gene_metrics %>%
#   ggplot(aes(x = rank_invar, y = rank_invar_velm, color = Gene.Type)) +
#   geom_text(aes(label = ifelse(!is.na(Gene.Type),Symbol,"")), color = "black") +
#   geom_point()+
#   theme_bw()

scatter_RI <- ggpairs(gene_metrics,
    columns = c("rank_invar", "rank_invar_velm", "rank_invar_velm_broad", "rank_invar_velm_age")
)

ggsave(scatter_RI, filename = here(plot_dir, "explore", "velm_RI_scatter.png"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
