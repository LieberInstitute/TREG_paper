library("SingleCellExperiment")
library("tidyverse")
library("GGally")
library("here")
library("TREG")

#### Load Data ####
gene_metrics <- read.csv(here("processed-data", "01_find_tregs", "supp_tables", "gene_metrics2.csv"), row.names = 1)

## Mathys data
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/referenceDatasets/SCE_mathys-PFC-BA10_MNT.rda", verbose = TRUE)
dim(sce.mathys)
# [1] 17923 70634

table(sce.mathys$Dx)
# Control      AD
# 44486   26148


## Populate rowData
## Only have Symbols
rowData(sce.mathys)

table(rownames(sce.mathys) %in% gene_metrics$Symbol)
# FALSE  TRUE
# 2520 15403

#### Filter to Prop Zero Genes from 10x data ####
gene_set <- gene_metrics$Symbol[gene_metrics$PropZero_filter]
length(gene_set)
# [1] 877
table(gene_set %in% row.names(sce.mathys))
# FALSE  TRUE
#    37   840
sce.mathys <- sce.mathys[gene_set[gene_set %in% row.names(sce.mathys)], ]
dim(sce.mathys)


## Check out cell types
table(sce.mathys$broad.cell.type)
# Ast   End    Ex    In   Mic   Oli   Opc   Per
# 3392   121 34976  9196  1920 18235  2627   167
#
summary(sce.mathys$age_death)

## Run RI
rank_invar_mathys <- rank_invariance_express(sce.mathys, group_col = "broad.cell.type")
gene_metrics_mathys <- as.data.frame(rank_invar_mathys) %>%
    rownames_to_column("Symbol")

## Only use controls
sce.mathys <- sce.mathys[, sce.mathys$Dx == "Control"]
dim(sce.mathys)
# [1] 17923 44486

rank_invar_mathys_filter <- rank_invariance_express(sce.mathys, group_col = "broad.cell.type")
gene_metrics_mathys_filter <- as.data.frame(rank_invar_mathys_filter) %>%
    rownames_to_column("Symbol")

## Add to gene  metrics
gene_metrics <- gene_metrics %>%
    left_join(gene_metrics_mathys) %>%
    left_join(gene_metrics_mathys_filter)

gene_metrics %>%
    filter(rank_invar > 850) %>%
    arrange(-rank_invar)

gene_metrics %>%
    filter(rank_invar < 10) %>%
    arrange(-rank_invar)

write.csv(gene_metrics, file = here("processed-data", "01_find_tregs", "gene_metrics_mathys.csv"))

#### Plots ####
plot_dir <- "plots/01_find_tregs"

# scatter_RI <- gene_metrics %>%
#   ggplot(aes(x = rank_invar, y = rank_invar_mathys, color = Gene.Type)) +
#   geom_text(aes(label = ifelse(!is.na(Gene.Type),Symbol,"")), color = "black") +
#   geom_point()+
#   theme_bw()

scatter_RI <- ggpairs(gene_metrics,
    columns = c("rank_invar", "rank_invar_mathys", "rank_invar_mathys_filter")
)

ggsave(scatter_RI, filename = here(plot_dir, "explore", "mathys_RI_scatter.png"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
