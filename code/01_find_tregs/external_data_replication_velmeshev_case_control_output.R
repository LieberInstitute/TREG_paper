
library("tidyverse")
library("here")
library("UpSetR")
library("org.Hs.eg.db")

## set up

plot_dir <- here("plots","01_find_tregs","external_data_replication_velmeshev_case_control")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

#### Load Data ####
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
  dplyr::slice(1)

nrow(rd) == nrow(sce.asd)

rowData(sce.asd) <- DataFrame(rd)
rowData(sce.asd)

## RI output
load(here("processed-data","01_find_tregs","velm_case_control_rank_invar.Rdata"), verbose = TRUE)
# velm_case_control_rank_invar

names(velm_case_control_rank_invar)

velm_case_control_rank_invar <- map(velm_case_control_rank_invar,
                                    ~.x |> left_join(rd) |>
                                      mutate(rank = row_number()))

map(velm_case_control_rank_invar,head)

top50 <- map(velm_case_control_rank_invar, ~head(.x, 50) |> pull(SYMBOL))

top50

pdf(file = here(plot_dir, "Velm_case_control_top50_TREG_upset.pdf"))
upset(fromList(top50), order.by = "freq")
dev.off()

## in velm data
Reduce("intersect", top50)
# [1] "MALAT1"    "LINC00486" "RASGEF1B"  "GPM6A"     "NCOR1"     "CDK13"     "ROBO2"     NA          "IL1RAPL1" 
# [10] "LRP1B"     "CADM2"     "LSAMP"     "ADGRB3"    "DLGAP1"    "NRXN1"     "HERC4"     "CHCHD3"    "SENP6"    
# [19] "BRWD1"     "SPPL3"     "NBEA"      "PCDH9"     "MEG3"      "LRRTM4"    "VTI1A"     "ARID4B"    "DENND4A"  
# [28] "COP1"      "RBM26"     "PUM2"      "TMCC1"     "DLG2"      "SRSF11"    "NRG3"     

map(top50, ~"AR1D1B" %in% .x)

map_dfr(velm_case_control_rank_invar,~.x |> filter(SYMBOL == "CADM2")) |>
  add_column(set = names(velm_case_control_rank_invar), .before = 1)

#       set         ENSEMBL rank_invar_velm SYMBOL rank
# 1     all ENSG00000175161             424  CADM2   11
# 2    case ENSG00000175161             391  CADM2   14
# 3 control ENSG00000175161             463  CADM2    7

## AKt3 & ARID1B
map_dfr(velm_case_control_rank_invar,~.x |> filter(SYMBOL %in% c("AKT3","ARID1B")))
#           ENSEMBL rank_invar_velm SYMBOL rank
# 1 ENSG00000117020              75   AKT3  360
# 2 ENSG00000049618              27 ARID1B  408
# 3 ENSG00000117020              78   AKT3  327
# 4 ENSG00000049618              22 ARID1B  383
# 5 ENSG00000117020              72   AKT3  398
# 6 ENSG00000049618              26 ARID1B  444

library(ggvenn)

## what about our TREG data
load(here("processed-data", "01_find_tregs", "rank_invar.Rdata"), verbose = TRUE)
# gene_propZero
# rank_invar_df

treg_top50 <- head(rank_invar_df, 50) |>
  mutate(rank = row_number())

Reduce("intersect", c(top50, list(treg = unlist(treg_top50$Symbol))))
# [1] "MALAT1" "CADM2"


pdf(file = here(plot_dir, "Velm_case_control_top50_TREG_Venn.pdf"), height = 5, width = 5)
ggvenn(
  top50, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()

#### Plot case v control ####
genes_of_interest <- rd |> filter(SYMBOL %in% c("MALAT1", "CADM2")) 

# hk_counts <- log2(as.matrix(assays(sce.asd[genes_of_interest$ENSEMBL, ])$counts) + 1)
## no counts :(
hk_counts <- assays(sce.asd[genes_of_interest$ENSEMBL, ])$logcounts

pd <- colData(sce.asd) |>
  as.data.frame() |>
  dplyr::select(ID = cell, cluster, sample, sum = UMIs, diagnosis )

counts_long <- as.data.frame(hk_counts) |>
  rownames_to_column("ENSEMBL") |>
  pivot_longer(!ENSEMBL, names_to = "ID", values_to = "logcount") |>
  left_join(pd) |>
  mutate(logsum = log2(as.double(sum))) |>
  left_join(genes_of_interest)


## plot
hk_sum_smooth <- counts_long |>
  filter(SYMBOL == "CADM2") |>
  ggplot() +
  geom_smooth(aes(logsum, logcount, color = cluster), method = "lm") +
  facet_wrap(~diagnosis) +
  theme_bw() +
  labs(x = "log2(sum)")

ggsave(hk_sum_smooth, filename = here(plot_dir,  "Velm_case_control_hk_sum_smooth.png"), width = 10, height = 5)

hk_sum_scater <- counts_long |>
  filter(SYMBOL == "CADM2") |>
  ggplot(aes(logsum, logcount, color = cluster)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~diagnosis) +
  theme_bw() +
  labs(x = "log2(sum)")

ggsave(hk_sum_scater, filename = here(plot_dir,  "Velm_case_control_hk_sum_scater.png"), width = 10, height = 5)



## linear modeling ####
mod_ct <- model.matrix(~ log2(sum) + cellType.Broad, data = colData(sce.asd))

fit <- lmFit(log2(assays(sce_pan)$counts + 1), mod_ct)
eB <- eBayes(fit)
tt <- topTable(eB, coef = "log2(sum)", number = Inf, sort.by = "none")

