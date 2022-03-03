library(SingleCellExperiment)
library(jaffelab)
library(edgeR)
library(limma)
library(here)

#### Prep sce data ####
load(here("raw-data", "sce_pan.v2.Rdata"), verbose = TRUE)
summary(sce_pan$sum)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 102    5766    9813   17282   23563  196431 

#### Find Expression Rank ####
## takes like 50 min
rank_df <- apply(as.matrix(assays(sce_pan)$logcounts), 2, rank) %>%
  as.data.frame()

save(rank_df, file = here("processed-data", "01_find_tregs", "rank_df.Rdata"))

## subset to canididate genes "AKT3"   "ARID1B" "MALAT1" "POLR2A"
rank_df_subset <- rank_df[c("ENSG00000117020", "ENSG00000049618", "ENSG00000251562", "ENSG00000181222"),]
save(rank_df_subset, file = here("processed-data", "01_find_tregs", "rank_df_subset.Rdata"))

#### Linear modeling ####
mod_ct <- model.matrix(~log2(sum) + cellType.Broad, data = colData(sce_pan))

fit = lmFit(log2(assays(sce_pan)$counts + 1), mod_ct)
eB = eBayes(fit)
tt = topTable(eB, coef= "log2(sum)", number= Inf, sort.by = "none")
table(tt$adj.P.Val < 0.05)
# FALSE  TRUE 
#    14 23024 

save(tt, file = here("processed-data", "01_find_tregs", "lmfit.Rdata"))

# sgejobs::job_single('snRNA_validation', create_shell = TRUE, memory = '100G', command = "Rscript snRNA_validation.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
