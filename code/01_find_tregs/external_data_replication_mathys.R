library("SingleCellExperiment")
library("tidyverse")
library("here")
library("TREG")

#### Load Data ####
gene_metrics <- read.csv(here("processed-data", "01_find_tregs", "supp_tables", "gene_metrics2.csv"), row.names = 1)

## Mathys data
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/referenceDatasets/SCE_mathys-PFC-BA10_MNT.rda", verbose = TRUE)
dim(sce.mathys)
# [1] 17923 70634

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
sce.mathys <- sce.mathys[gene_set[gene_set %in% row.names(sce.mathys)],]
dim(sce.mathys)

# ## add ensmbl ids 
# rd <- rowData(sce.mathys) %>% 
#   as.data.frame() %>% 
#   rename(Symbol = V1) %>%
#   left_join(gene_metrics %>% select(Symbol, ensembl_id, PropZero_filter))
# 
# rowData(sce.mathys) <- DataFrame(rd)
# 
# ## Check out cell types
# table(sce.mathys$broad.cell.type)
# # Ast   End    Ex    In   Mic   Oli   Opc   Per 
# # 3392   121 34976  9196  1920 18235  2627   167
# 
# # sce.mathys <- sce.mathys[,!sce.mathys$broad.cell.type %in% c("End", "Per")]
# dim(sce.mathys)
# # [1]   840 70346

## Run RI
rank_invar_mathys <- rank_invariance_express(sce.mathys, group_col = "broad.cell.type")
gene_metrics_mathys <- as.data.frame(rank_invar_mathys) %>%
  rownames_to_column("Symbol")
  
## Add to gene  metrics
gene_metrics <- gene_metrics %>%
  left_join(gene_metrics_mathys)

gene_metrics %>%
  filter(rank_invar > 850) %>%
  arrange(-rank_invar)

gene_metrics %>%
  filter(rank_invar <10) %>%
  arrange(-rank_invar) 

#### Plots ####
plot_dir <- "plots/01_find_tregs"

scatter_RI <- gene_metrics %>% 
  ggplot(aes(x = rank_invar, y = rank_invar_mathys, color = Gene.Type)) +
  geom_text(aes(label = ifelse(!is.na(Gene.Type),Symbol,"")), color = "black") +
  geom_point()+
  theme_bw()

ggsave(scatter_RI, filename = here(plot_dir, "explore","mathys_RI_scatter.png"))

scatter_diff_RI <- gene_metrics %>% 
  ggplot(aes(x = rank_invar, y = rank_invar - rank_invar_mathys, color = Gene.Type)) +
  geom_point() 

ggsave(scatter_diff_RI, filename = here(plot_dir, "explore","mathys_RI_diff_scatter.png"))

density_diff_RI <- gene_metrics %>% 
  ggplot(aes(x = rank_invar - rank_invar_mathys)) +
  geom_density() 

ggsave(density_diff_RI, filename = here(plot_dir, "explore","mathys_RI_diff_density.png"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
