library("tidyverse")
library("here")
library("sessioninfo")
library("RColorBrewer")
library("ggrepel")

gene_metrics <- read.csv(here("processed-data", "01_find_tregs","supp_tables", "gene_metrics2.csv"), row.names = 1)
head(gene_metrics)

## Marker genes from TRan Maynard et al. 
## https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/810b47364af4c8afe426bd2a6b559bd6a9f1cc98/tables/revision/top40genesLists_DLPFC-n3_cellType_SN-LEVEL-tests_LAH2020.csv
tran_markers <- read.csv(here("processed-data","00_data_prep","top40genesLists_DLPFC-n3_cellType_SN-LEVEL-tests_LAH2020.csv"))

tran_markers_long <- tran_markers %>%
  rownames_to_column("rank") %>%
  pivot_longer(cols = !rank, names_to = "Marker", values_to = "Symbol") %>%
  filter(Symbol != "") %>%
  separate(Marker, into = c("cellType", "test"), "_", extra = "merge", remove = FALSE) %>%
  separate(test, into = c("cellType_suf", "test"), "_", fill = "left") %>%
  mutate(cellType = ifelse(is.na(cellType_suf), cellType, paste0(cellType, "_" ,cellType_suf)),
         test = paste0("Tran_", test)) %>%
  select(-cellType_suf) 

tran_markers_long %>% count(cellType, test)

# tran_markers_wide <- tran_markers_long %>%
#   pivot_wider(names_from = "test", values_from = tran_marker, names_prefix = "marker_tran_") %>%
#   replace_na(list(marker_tran_1vAll = FALSE, marker_tran_pw = FALSE))
# 

## Make table for Park markers
park_markers_list <- readLines(here("raw-data","Park_markers_n69.txt"))
park_markers <- tibble(Symbol = park_markers_list, 
                       test = "Park")

markers_anno <- tran_markers_long %>%
  select(Symbol, test) %>%
  rbind(park_markers) %>%
  group_by(Symbol) %>%
  mutate(test = ordered(test, levels = c("Tran_1vAll","Tran_pw","Park"))) %>%
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

gene_metrics_marker <- gene_metrics %>% 
  left_join(markers_anno, by = "Symbol") %>%
  mutate(`Gene Type` = factor(ifelse(!is.na(Marker), Marker, Gene.Type),
                              levels = c("TREG Canidate", "Classic HK", "Data Driven HK",
                                         "Tran_1vAll", "Tran_pw", "Park", "Tran_1vAll,Tran_pw",
                                         "Tran_1vAll,Park", "Tran_1vAll,Tran_pw,Park", "None"))) %>%
  replace_na(list(`Gene Type` = "None"))

## no overlapping Marker + Gene type annotations
gene_metrics_marker %>%
  count(Marker, Gene.Type, `Gene Type`) %>%
  arrange(`Gene Type`)

## 11 genes are markers + have RI vals
gene_metrics_marker %>%
  filter(!is.na(Marker) & !is.na(rank_invar)) %>%
  arrange(-rank_invar)

#     Symbol      ensembl_id top50 max_PropZero PropZero_filter rank_invar Gene.Type         t             Marker
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


#### Linear modeling ####
# load(here("processed-data", "01_find_tregs", "lmfit.Rdata"), verbose = TRUE) # tt

## bind with rank invar data
gene_metrics_marker2 <- gene_metrics_marker %>%
  mutate(
    # Signif = cut(adj.P.Val, c(0, 0.001 ,0.01, 0.05, 1), include.lowest = TRUE),
         gene_anno = case_when(
           !top50 ~ "Fail 50% Exp.",
           !PropZero_filter ~ "Fail Prop. Zero",
           TRUE ~ "RI"),
         label1 = Symbol %in% c("ARID1B","AKT3","MALAT1","POLR2A"),
         label2 = `Gene Type` != "None",
         alpha = `Gene Type` != "None"
  )

gene_metrics_marker2 %>% dplyr::count(gene_anno, `Gene Type`)
# gene_anno               Gene Type     n
# 1    Fail 50% Exp.              Classic HK     2
# 2    Fail 50% Exp.          Data Driven HK     1
# 3    Fail 50% Exp.              Tran_1vAll    14
# 4    Fail 50% Exp.                 Tran_pw    24
# 5    Fail 50% Exp.                    Park    14
# 6    Fail 50% Exp.      Tran_1vAll,Tran_pw   100
# 7    Fail 50% Exp.         Tran_1vAll,Park     1
# 8    Fail 50% Exp. Tran_1vAll,Tran_pw,Park     3
# 9    Fail 50% Exp.                    None 11360
# 10 Fail Prop. Zero              Classic HK    10
# 11 Fail Prop. Zero          Data Driven HK     7
# 12 Fail Prop. Zero              Tran_1vAll   186
# 13 Fail Prop. Zero                 Tran_pw   158
# 14 Fail Prop. Zero                    Park    17
# 15 Fail Prop. Zero      Tran_1vAll,Tran_pw   292
# 16 Fail Prop. Zero         Tran_1vAll,Park     7
# 17 Fail Prop. Zero Tran_1vAll,Tran_pw,Park    17
# 18 Fail Prop. Zero                    None  9948
# 19              RI           TREG Canidate     3
# 20              RI              Tran_1vAll     3
# 21              RI                 Tran_pw     6
# 22              RI      Tran_1vAll,Tran_pw     2
# 23              RI                    None   863

## t-stat plots
plot_dir <- "plots/01_find_tregs"
gene_type_colors <- c(brewer.pal(n = 9, name = "Spectral"), "grey")
# gene_type_colors <- c(brewer.pal(n = 9, name = "Spectral"), "black")
names(gene_type_colors) <- levels(gene_metrics_marker$`Gene Type`)

invar_t_scatter <- gene_metrics_marker2 %>%
  filter(gene_anno == "RI") %>% 
  ggplot(aes(x = rank_invar, y = t))+
  geom_point(alpha = 0.2, color = "grey") +
  geom_point(data = filter(gene_metrics_marker2, label2, gene_anno == "RI"), 
             aes(fill= `Gene Type`), shape=21, color = "black") +
  geom_text_repel(aes(label = ifelse(label2, paste0("italic('",Symbol,"')"), NA)),
                  size = 3, parse = TRUE) +
  scale_fill_manual(values = gene_type_colors, drop = FALSE) +
  ylim(0,400) +
  theme_bw() +
  labs(x = "Rank Invariance") +
  theme(text = element_text(size=15))

ggsave(invar_t_scatter, filename = here(plot_dir, "explore" ,"rank_invar_t_scatter.png"), width = 6, height = 6.3)
ggsave(invar_t_scatter, filename = here(plot_dir, "supp_pdf" ,"rank_invar_t_scatter.pdf"), width = 6, height = 6.3)

## Denisity/jitter plot
pos <- position_jitter(width = 0.3, seed = 2)
invar_t_density <- gene_metrics_marker2 %>%
  filter(gene_anno != "RI") %>%
  ggplot(aes(x = `Gene Type`, y = t))+
  geom_point(aes(color = `Gene Type`, alpha = alpha), position = pos) +
  # geom_point(aes(alpha = alpha, color = `Gene Type`), position = pos) +
  # geom_point(data = filter(gene_metrics_marker2, label2, gene_anno != "RI"), 
  #            aes(fill= `Gene Type`), shape=21, color = "black", position = pos)+
  geom_text_repel(aes(label = ifelse(label1, paste0("italic('",Symbol,"')"), NA)), 
                  size = 3, position = pos, parse = TRUE) +
  scale_color_manual(values = gene_type_colors, drop = FALSE) +
  facet_wrap(~gene_anno)+
  ylim(0,400) +
  theme_bw()+ 
  theme(legend.position="none",
        text = element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x=element_blank()
  ) 

ggsave(invar_t_density , filename = here(plot_dir,"explore", "rank_invar_t_densitiy.png"), width = 5)
ggsave(invar_t_density , filename = here(plot_dir,"supp_pdf", "rank_invar_t_densitiy.pdf"), width = 5)

## add to gene metrics and export
gene_metrics2 <- gene_metrics %>%
  left_join(tt %>% select(t) %>% rownames_to_column("ensembl_id"))

write.csv(gene_metrics2, here("processed-data", "01_find_tregs","supp_tables", "gene_metrics2.csv"))

