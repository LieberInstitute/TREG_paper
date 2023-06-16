library("SingleCellExperiment")
library("tidyverse")
library("here")
# library("TREG")

## get Dx input
args <- commandArgs(trailingOnly = TRUE)
dx_i <- as.integer(args[[1]])

## Dx sets

dx_sets <- list(all = c("ASD", "Control"), case = "ASD", control = "Control")
dx_name <- names(dx_sets)[[dx_i]]
message(Sys.time(), " - Load data, will subset to: ", dx_name)

## Load Data
cell_type_match <- read.csv(here("processed-data", "01_find_tregs", "Velmeshev_broad.csv")) |>
    dplyr::select(-Note)

load("/dcl01/ajaffe/data/lab/singleCell/velmeshev2019/analysis_MNT/SCE_asd-velmeshev-etal_MNT.rda",
     verbose = T
) # sce.asd

#### Match Up Cell Types ####
pd <- colData(sce.asd) |>
    as.data.frame() |>
    left_join(cell_type_match)

colData(sce.asd) <- DataFrame(pd)

## no cell types < 50 
table(sce.asd$cellType.Broad, sce.asd$diagnosis)
# Astro  Endo Excit Inhib Micro Oligo   OPC 
# 12217  2668 44355 16494  3331 15371 10123


## rank
sce.asd <- sce.asd[, sce.asd$diagnosis %in% dx_sets[[dx_i]]]
message("subset to ", dx_name, ", ncol  = ", ncol(sce.asd))

message(Sys.time(), " - Ranking genes")

rank_df <- apply(as.matrix(assays(sce.asd)$logcounts), 2, rank) %>%
  as.data.frame()

message(Sys.time(), " - save data")
save(rank_df, file = here("processed-data", "01_find_tregs", paste0("velm_rank_df_",dx_name,".Rdata")))

# sgejobs::job_single('external_data_replication_velmeshev_get_rank', create_shell = TRUE, memory = '100G', command = "Rscript external_data_replication_velmeshev_get_rank.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
