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
eb_HKG <- c("C1orf43", "CHMP2A", "EMC7", "GPI", "PSMB2","PSMB4","RAB7A","REEP5","SNRPD3","VCP","VPS29")
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


gene_metrics %>%filter(Symbol %in% eb_HKG)

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
        `Gene Type` = factor(`Gene Type`, levels = c('TREG Candidate', 'Classic HK', "Data Driven HK", "None"))
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
    geom_text_repel(aes(label = ifelse(label2, paste0("italic('", Symbol, "')"), NA),
                        color = label1),
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
    geom_text_repel(aes(label = ifelse(label2, paste0("italic('", Symbol, "')"), NA),
                        color = label1),
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
all_entrez <- bitr(gene_metrics2$ensembl_id, fromType = "ENSEMBL", toType ="ENTREZID", OrgDb="org.Hs.eg.db")
nrow(all_entrez)
# [1] 18296
u <- all_entrez$ENTREZ

## List top 10% (87) of the rank invar genes
top_RI <- gene_metrics2 %>%
    filter(!is.na(rank_invar)) %>%
    arrange(-rank_invar) %>%
    head(sum(!is.na(gene_metrics2$rank_invar))%/%10) %>%
    pull(ensembl_id)

top_RI_entrez <- bitr(top_RI, fromType = "ENSEMBL", toType ="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZ
length(top_RI_entrez)
# [1] 86

## Run Enrichment 

# go_all <- compareCluster(geneClusters = list(topRI = top_RI_entrez), 
#                          univ = u,
#                          OrgDb = "org.Hs.eg.db", 
#                          fun = "enrichGO",
#                          ont = "ALL")


ont <- c("BP", "CC", "MF", "ALL")
names(ont) <- ont
go_output = map(ont, ~compareCluster(geneClusters = list(topRI = top_RI_entrez),
                                       univ = u,
                                       OrgDb = "org.Hs.eg.db",
                                       fun = "enrichGO",
                                       ont = .x))

pdf(here(plot_dir, "explore", "go_all_top87.pdf"))
map2(go_output, names(go_output), ~dotplot(.x, title = .y))
dev.off()

save(go_output, file = here("processed-data", "01_find_tregs","enrichGO.Rdata"))

map(go_output, as.data.frame)
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
    left_join(marker_stats %>% select(Symbol, cellType))  %>%
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
