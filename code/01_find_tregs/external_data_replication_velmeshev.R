library("SingleCellExperiment")
# library(scater)
library("tidyverse")
library("here")
library("TREG")
## Load Data
gene_metrics <- read.csv(here("processed-data", "01_find_tregs", "supp_tables", "gene_metrics2.csv"), row.names = 1)

load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/SCE_asd-velmeshev-etal_controls_multiBatchNorm.Rdata", verbose = TRUE)
dim(sce.asd)
# 36501 52556

## Populate rowData
rowData(sce.asd)

(rowData(sce.asd) <- DataFrame(ensmblID = rownames(sce.asd),
                               Symbol = guess_symbol))
#### Filter to Prop Zero Genes from 10x data ####
gene_set <- gene_metrics$ensembl_id[gene_metrics$PropZero_filter]
length(gene_set)
# [1] 877
table(gene_set %in% row.names(sce.asd))
# FALSE  TRUE 
#     1   876
sce.asd <- sce.asd[gene_set[gene_set %in% row.names(sce.asd)],]
dim(sce.asd)

table(sce.asd$cluster)
levels(sce.asd$cluster)

## Run RI
rank_invar_velm <- rank_invariance_express(sce.asd, group_col = "cluster")
gene_metrics_velm <- as.data.frame(rank_invar_velm) %>%
  rownames_to_column("ensembl_id")
  

## Add to gene  metrics
gene_metrics <- gene_metrics %>%
  left_join(gene_metrics_velm)


gene_metrics %>%
  filter(rank_invar > 867) %>%
  arrange(-rank_invar)

gene_metrics %>%
  filter(rank_invar <10) %>%
  arrange(-rank_invar) 

#### Plots ####
plot_dir <- "plots/01_find_tregs"

scatter_RI <- gene_metrics %>% 
  ggplot(aes(x = rank_invar, y = rank_invar_velm, color = Gene.Type)) +
  geom_text(aes(label = ifelse(!is.na(Gene.Type),Symbol,"")), color = "black") +
  geom_point()+
  theme_bw()

ggsave(scatter_RI, filename = here(plot_dir, "explore","velm_RI_scatter.png"))

scatter_diff_RI <- gene_metrics %>% 
  ggplot(aes(x = rank_invar, y = rank_invar - rank_invar_velm, color = Gene.Type)) +
  geom_point() 

ggsave(scatter_diff_RI, filename = here(plot_dir, "explore","velm_RI_diff_scatter.png"))

density_diff_RI <- gene_metrics %>% 
  ggplot(aes(x = rank_invar - rank_invar_velm)) +
  geom_density() 

ggsave(density_diff_RI, filename = here(plot_dir, "explore","velm_RI_diff_density.png"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
