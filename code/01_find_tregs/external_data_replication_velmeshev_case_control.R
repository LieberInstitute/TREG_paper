library("SingleCellExperiment")
library("tidyverse")
library("here")
library("TREG")
library("org.Hs.eg.db")

plot_dir <- here("plots", "01_find_tregs", "external_data_replication")
if(!dir.exists(plot_dir)) dir.create(plot_dir)
  
## Load Data
# gene_metrics <- read.csv(here("processed-data", "01_find_tregs", "supp_tables", "gene_metrics2.csv"), row.names = 1)
cell_type_match <- read.csv(here("processed-data", "01_find_tregs", "Velmeshev_broad.csv")) %>%
    dplyr::select(-Note)

load("/dcl01/ajaffe/data/lab/singleCell/velmeshev2019/analysis_MNT/SCE_asd-velmeshev-etal_MNT.rda",
     verbose = T
) # sce.asd
sce.asd
# 36501 104559

table(sce.asd$diagnosis, sce.asd$region)
#           ACC   PFC
# ASD     19984 32019
# Control 22409 30147

## Populate rowData
annots <- select(org.Hs.eg.db, keys=rownames(sce.asd), 
                 columns="SYMBOL", keytype="ENSEMBL")
nrow(annots)
# 'select()' returned 1:many mapping between keys and columns
annots |> group_by(ENSEMBL) |> filter(n()  > 2)

## just pick first match - not perfect...
rd <- annots |>
  group_by(ENSEMBL) |>
  slice(1)

nrow(rd) == nrow(sce.asd)

rowData(sce.asd) <- DataFrame(rd)
rowData(sce.asd)

rd |> filter(SYMBOL %in% c("MALAT1", "AKT3", "ARID1B"))

## note most donors are < 17
table(sce.asd$age)
# 4     6     7     8    11    12    13    14    15    16    19    20    21    22 
# 2973  2255  5915  2176  2812 17343 11528 10033 14720  1112 15467  4670  8871  4684

#### Match Up Cell Types ####
pd <- colData(sce.asd) %>%
    as.data.frame() %>%
    left_join(cell_type_match)

colData(sce.asd) <- DataFrame(pd)

unique(sce.asd$cluster)
table(sce.asd$cluster, sce.asd$cellType.Broad)
#                  Astro  Endo Excit Inhib Micro Oligo   OPC
# AST-FB            4562     0     0     0     0     0     0
# AST-PP            7655     0     0     0     0     0     0
# Endothelial          0  2668     0     0     0     0     0
# IN-PV                0     0     0  4037     0     0     0
# IN-SST               0     0     0  4569     0     0     0
# IN-SV2C              0     0     0  1914     0     0     0
# IN-VIP               0     0     0  5974     0     0     0
# L2/3                 0     0 13324     0     0     0     0
# L4                   0     0  6858     0     0     0     0
# L5/6                 0     0  3561     0     0     0     0
# L5/6-CC              0     0  4617     0     0     0     0
# Microglia            0     0     0     0  3331     0     0
# Neu-mat              0     0  4135     0     0     0     0
# Neu-NRGN-I           0     0  3543     0     0     0     0
# Neu-NRGN-II          0     0  8317     0     0     0     0
# Oligodendrocytes     0     0     0     0     0 15371     0
# OPC                  0     0     0     0     0     0 10123

## no cell types < 50 
table(sce.asd$cellType.Broad, sce.asd$diagnosis)
# Astro  Endo Excit Inhib Micro Oligo   OPC 
# 12217  2668 44355 16494  3331 15371 10123

## Dx sets

dx_sets <- list(all = c("ASD", "Control"), case = "ASD", control = "Control")

velm_case_control_rank_invar <- map2(dx_sets, names(dx_sets), function(dx, name){
  
  ## subset 
  
  sce.asd <- sce.asd[, sce.asd$diagnosis %in% dx]
  message("subset to ", name, ", nrow  = ", nrow(sce.asd))
  
  #### Filter to top 50% and low Prop Zero Genes from 10x data ####
  
  ## filter to top 50%
  message(Sys.time(), " - Top 50% filter")
  row_means <- rowMeans(assays(sce.asd)$logcounts)
  (median_row_means <- median(row_means))
  # [1] 0.0121892
  
  sce.asd <- sce.asd[row_means > median_row_means, ]
  message("Subset to top 50%: nrow:", nrow(sce.asd))
  # [1] 18250
  
  ## prop zero
  ## can't use logcounts in get_prop_zero fake it - fix this in package
  names(assays(sce.asd)) <- "counts"
  
  message(Sys.time(), " - Get Prop Zero")
  prop_zeros <- get_prop_zero(sce.asd, group_col = "cellType.Broad")
  head(prop_zeros)
  
  ## swap it back
  names(assays(sce.asd)) <- "logcounts"
  
  ## plot prop zero
  # Pivot data longer for plotting
  prop_zero_long <- prop_zeros %>%
    rownames_to_column("Gene") %>%
    pivot_longer(!Gene, names_to = "Group", values_to = "prop_zero")
  
  # Plot histograms
  propZero_limit = 0.85
  
  prop_zero_histogram <- ggplot(
    data = prop_zero_long,
    aes(x = prop_zero, fill = Group)
  ) +
    geom_histogram(binwidth = 0.05) +
    facet_wrap(~Group) +
    geom_vline(xintercept = propZero_limit, color = "red", linetype = "dashed")
  
  ggsave(prop_zero_histogram, filename = here(plot_dir,paste0("Velm_",name,"_prop_zero_histo.png")))
  
  ## filter genes 
  message(Sys.time(), "filter prop zero")
  filtered_genes <- filter_prop_zero(prop_zeros, cutoff = propZero_limit)
  message("n genes: ",length(filtered_genes))
  # [1] 434 with limit = 0.85
  
  ## check 10x TREGs 
  # rd |> filter(SYMBOL %in% c("MALAT1", "AKT3", "ARID1B")) |> mutate(ENSEMBL %in% filtered_genes)
  
  ## filter
  # length(filtered_genes) / nrow(sce.asd)
  sce.asd <- sce.asd[filtered_genes, ]
  
  
  ## Run RI with all data, broad cell types
  rank_invar_velm <- rank_invariance_express(sce.asd, group_col = "cellType.Broad")
  
  gene_metrics_velm_broad <- as.data.frame(rank_invar_velm) |>
    rownames_to_column("ENSEMBL") |> 
    arrange(-rank_invar_velm_all_broad) |>
    left_join(rd)
  
  # head(gene_metrics_velm_broad, n = 25)
  
  # gene_metrics_velm_broad |> filter(SYMBOL %in% c("MALAT1", "AKT3", "ARID1B"))
  
  return(rank_invar_velm)
  
})

save(velm_case_control_rank_invar, here("processed_data","01_find_tregs","velm_case_control_rank_invar.Rdata"))

# sgejobs::job_single('external_data_replication_velmeshev_case_control', create_shell = TRUE, memory = '100G', command = "Rscript external_data_replication_velmeshev_case_control.R")


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
