library(SingleCellExperiment)
library(SummarizedExperiment)
library(EvenExpress)
library(tidyverse)
library(here)

load(here("data","sce_pan.v2.Rdata"), verbose = TRUE)

table(sce_pan$region, sce_pan$cellType.Broad)
# orignal_regions <- sce_pan$region ## prehaps look at rare cell types in dlpfc?

sce_pan$region[sce_pan$cellType.Broad %in% c("Macro", "Mural", "Endo", "Tcell")] <- "combined"
sce_pan$ctXregion <-paste0(sce_pan$cellType.Broad, "_" ,sce_pan$region)

## Calc sum_counts
sce_pan$sum_counts <- colSums(assays(sce_pan)$counts)

## filter for top 50% expressed of genes
row_means <- rowMeans(assays(sce_pan)$logcounts)
median_row_means <- median(row_means)

sce_pan <- sce_pan[row_means > median_row_means,]
dim(sce_pan)

gene_propZero <- get_prop_zero(sce_pan, "ctXregion")
head(gene_propZero)

propZero_limit <- 0.75
# genes_filtered <- filter_prop_zero(gene_propZero, cutoff = propZero_limit)
max_prop_zero <- apply(gene_propZero, 1, max)
genes_filtered <- names(max_prop_zero[max_prop_zero < propZero_limit])
length(genes_filtered)
# [1] 877

sce_pan <- sce_pan[genes_filtered,]
dim(sce_pan)
# [1]   877 70527

## Run rank invar for full data set
assays(sce_pan)$logcounts <- as.matrix(assays(sce_pan)$logcounts)
gene_symbol <- as.data.frame(rowData(sce_pan)) %>% select(gene = gene_id, Symbol)

rank_invariance <- get_rank_invariance_express(sce_pan, "cellType.Broad")

rank_invar_df <- as.data.frame(rank_invariance)
colnames(rank_invar_df ) <- "rank_invar"
rank_invar_df <- rank_invar_df %>% rownames_to_column("gene") %>% left_join(gene_symbol) %>% arrange(-rank_invar)

head(rank_invar_df, 20)
#               gene rank_invar     Symbol
# 1  ENSG00000251562        877     MALAT1
# 2  ENSG00000171988        876     JMJD1C
# 3  ENSG00000230590        875        FTX
# 4  ENSG00000127603        874      MACF1
# 5  ENSG00000117020        873       AKT3
# 6  ENSG00000100354        872     TNRC6B
# 7  ENSG00000285106        871 AC016831.7
# 8  ENSG00000156639        870     ZFAND3
# 9  ENSG00000049618        869     ARID1B
# 10 ENSG00000175161        868      CADM2
# 11 ENSG00000142599        867       RERE
# 12 ENSG00000124788        866      ATXN1
# 13 ENSG00000120071        865     KANSL1
# 14 ENSG00000075151        864     EIF4G3
# 15 ENSG00000115524        863      SF3B1
# 16 ENSG00000181722        862     ZBTB20
# 17 ENSG00000132155        861       RAF1
# 18 ENSG00000125676        860      THOC2
# 19 ENSG00000055609        859      KMT2C
# 20 ENSG00000100393        858      EP300

## run for dlpfc
sce_dlpfc <- sce_pan[,sce_pan$region == "dlpfc"]
dim(sce_dlpfc)
# [1]   877 11165
sce_dlpfc$cellType.Broad <- droplevels(sce_dlpfc$cellType.Broad)
table(sce_dlpfc$cellType.Broad)

rank_invariance_dlpfc <- get_rank_invariance_express(sce_dlpfc, "cellType.Broad")

rank_invar_df_dlpfc <- as.data.frame(rank_invariance_dlpfc)
colnames(rank_invar_df_dlpfc) <- "rank_invar"
rank_invar_df_dlpfc  <- rank_invar_df_dlpfc %>% rownames_to_column("gene") %>% left_join(gene_symbol) %>% arrange(-rank_invar)


save(gene_propZero, rank_invar_df, rank_invar_df_dlpfc, file = here("data","rank_invar.Rdata"))

# sgejobs::job_single('find_housekeeper_genes_EvenExpress', create_shell = TRUE, queue= 'bluejay', memory = '100G', command = "Rscript find_housekeeper_genes_EvenExpress.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
