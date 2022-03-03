library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(jaffelab)
library(VennDiagram)
library(ggrepel)
library(here)

#### Fig 1. Demo Rank Distribution ####
set.seed(10)
demo_rank_data <- tibble(sample = 1:100,
                         `high\ninvariance` = round(rnorm(100, 75, 3), digits = 0),
                         `low\ninvariance` = sample(1:100, 100, replace = TRUE))  %>%
  pivot_longer(!sample, names_to = "Gene", values_to = "Expression Rank")

rank_demo_violin <- ggplot(demo_rank_data, aes(x = Gene, y = `Expression Rank`)) +
  geom_violin(fill = "light grey") +
  theme_bw() +
  theme(text = element_text(size=15))

# ggsave(rank_demo_violin, filename = here("plots", "01_find_tregs", "main_pdf","fig1_rank_violin_demo.png"), width = 3, height = 6)
ggsave(rank_demo_violin, filename = here("plots", "01_find_tregs", "main_pdf","fig1_rank_violin_demo.pdf"), width = 3, height = 3)

#### Prep sce data ####
load(here("raw-data", "sce_pan.v2.Rdata"), verbose = TRUE)
load(here("processed-data", "01_find_tregs", "rank_invar.Rdata"), verbose = TRUE)
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)

table(sce_pan$cellType.Broad)

## Capitalize Regions
sce_pan$region2 <- toupper(sce_pan$region)
sce_pan$region2[sce_pan$region2 == "SACC"] <- "sACC"
sce_pan$region2[sce_pan$region2 == "NAC"] <- "NAc"
## Rare cell types 
sce_pan$region2[sce_pan$cellType.Broad %in% c("Macro", "Mural", "Endo", "Tcell")] <- "combined"
table(sce_pan$region2)

sce_pan$ctXregion <-paste0(sce_pan$cellType.Broad, "_" ,sce_pan$region2)

pd <- as.data.frame(colData(sce_pan))

#### Annotate Genes ####
gene_metrics <- read.csv(here("processed-data", "01_find_tregs","supp_tables", "gene_metrics.csv"), row.names = 1)

gene_symbol <- rowData(sce_pan) %>% 
  as.data.frame() %>%
  select(gene = gene_id, Symbol)

### TREGs Chosen from top 10 RI gene for probe availability
treg_list <- c("MALAT1", "AKT3", "ARID1B")

# classic HK genes
dotdotdot_genes <- c("POLR2A","PPIB","UBC","HPRT1","ACTB","TUBB3","BIN1","LDHA",
                     "GAPDH","PGK1","BHLHE22","CPLX2")

data_driven_HK <- c("NDUFB4", "NDUFB1", "GSTO1", "AMZ2", "POLR2I", "NDUFA3", "RRAGA", "POMP")
  
(candidate_genes <- gene_symbol %>%
  filter(Symbol %in% c(treg_list,"POLR2A")))


genes_of_interest <- tibble(Symbol = c(treg_list, dotdotdot_genes, data_driven_HK),
                              gene_anno = c(rep("TREG Canidate", length(treg_list)),
                                            rep("Classic HK", length(dotdotdot_genes)),
                                            rep("Data Driven HK", length(data_driven_HK))
                                            )
                            ) %>%
  left_join(gene_symbol)

## annotate gene metrics
gene_metrics <- gene_metrics %>%
  left_join(genes_of_interest %>% select(Symbol, `Gene Type` = gene_anno)) 

gene_metrics %>%
  filter(!is.na(`Gene Type`))

gene_metrics %>%
  filter(Symbol %in% candidate_genes$Symbol)

#### Proportion zero plots ####
propZero_limit <- 0.75
## overall density

gene_propZero_long <- gene_propZero %>%
  rownames_to_column("gene") %>%
  pivot_longer(!gene, values_to = "propZero") %>%
  separate(name, into = c("cellType","region")) 

gene_propZero_long$cellType <- factor(gene_propZero_long$cellType, levels = levels(sce_pan$cellType.Broad))
gene_propZero_long$region <- factor(gene_propZero_long$region, levels = c("AMY", "DLPFC", "HPC", "NAc", "sACC", "combined"))

#### Fig 2 ####
propZero_density <- gene_propZero_long %>%
  filter(region != "combined") %>%
  ggplot(aes(x = propZero, fill = cellType)) +
  geom_histogram(binwidth = 0.05, color = "black",size=.2) +
  scale_fill_manual(values = cell_colors) +
  facet_grid(cellType~region) +
  labs(x = "Proportion Zero per Group", y = "Number of Genes")+
  geom_vline(xintercept = propZero_limit, color = "red", linetype = "dashed") +
  theme_bw()+
  theme(legend.position = "none", 
        text = element_text(size=15), 
        axis.text.x = element_text(angle = 90, hjust = 1))+ 
  scale_y_continuous(breaks=seq(0, 3000, 1500))+ 
  scale_x_continuous(breaks=seq(0, 1, .5))

# ggsave(propZero_density, filename = here("plots", "01_find_tregs", "main_pdf","fig2_propZero_density.png"), width = 6, height = 7)
ggsave(propZero_density, filename = here("plots", "01_find_tregs", "main_pdf","fig2_propZero_density.pdf"), width = 6, height = 7)


propZero_density_rare <- gene_propZero_long %>%
  filter(region == "combined") %>%
  ggplot(aes(x = propZero, fill = cellType)) +
  geom_histogram(binwidth = 0.05, color = "black",size=.2) +
  scale_fill_manual(values = cell_colors) +
  facet_wrap(~cellType, nrow = 1) +
  labs(x = "Proportion Zero per Group", y = "Number of Genes")+
  geom_vline(xintercept = propZero_limit, color = "red", linetype = "dashed") +
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=15), axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(propZero_density_rare, filename = here("plots", "01_find_tregs", "supp_pdf","fig_supp_propZero_density_rare_ct.png"), width = 7, height = 2)

## filter 
# sce_pan <- sce_pan[genes_filtered,]
# (n_genes <- c(n_genes, "zero_filter" = nrow(sce_pan)))
# sce_pan$sum_counts_zero_filter <- colSums(assays(sce_pan)$counts)
# dim(sce_pan)
# # [1]   877 70527

# propZero_density_filtered <- gene_propZero_long %>%
#   filter(gene %in% genes_filtered) %>%
#   ggplot(aes(x = propZero, fill = cellType)) +
#   geom_histogram(binwidth = 0.05, color = "black",size=.2) +
#   scale_fill_manual(values = cell_colors) +
#   facet_grid(region~cellType) +
#   geom_vline(xintercept = propZero_limit, color = "red", linetype = "dashed") +
#   theme_bw()+
#   theme(legend.position = "none", te) +
#   labs( title = paste0("Proportion Zero Distribution - filter for ", propZero_limit, " max non-zero"),
#         x = "Proportion Zero per Group", y = "Number of Genes")
# 
# ggsave(propZero_density_filtered, filename = here("plots", "01_find_tregs", "propZero_density_filter.png"), width = 10)

#### Demo filtering ####
filter_demo <- gene_symbol %>% 
  filter(Symbol %in% c(treg_list, "POLR2A")) %>%
  left_join(gene_propZero_long)

filter_anno <- filter_demo %>% 
  group_by(Symbol) %>%
  summarise(max = max(propZero))
  
filter_demo_scatter <- ggplot(filter_demo, aes(x = Symbol, y = propZero)) +
  geom_jitter(aes(color = cellType), width = 0.1) +
  scale_color_manual(values = cell_colors, name = "Cell Type") +
  labs(x = "Gene", y = "Proportion Zero") +
  geom_hline(yintercept = propZero_limit, color = "red", linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size=15), 
        axis.text.x = element_text(angle = 90, hjust = 1, face="italic"))

# ggsave(filter_demo_scatter, filename = here("plots", "01_find_tregs", "main_pdf","fig2b_propZero_filter_demo.png"), width = 4, height = 7)
ggsave(filter_demo_scatter, filename = here("plots", "01_find_tregs", "main_pdf","fig2b_propZero_filter_demo.pdf"), width = 4, height = 7)

#### Fig 3. Real Rank Distribution ####
load(here("processed-data", "01_find_tregs", "rank_df_subset.Rdata"), verbose = TRUE) 
# 4 canidiate genes * 70k nuclei
dim(rank_df_subset)
# [1]     4 70527 
corner(rank_df_subset)

rank_long <- rank_df_subset %>%
  rownames_to_column("gene") %>%
  pivot_longer(!gene, names_to = "Sample", values_to = "rank") %>%
  left_join(gene_symbol) %>%
  left_join(pd %>% select(Sample = uniqueID, cellType.Broad))

rank_violin <- ggplot(rank_long, aes(x = Symbol, y = rank)) +
  geom_violin(fill = "light grey", scale = "width") +
  labs(x = "Gene") +
  theme_bw() +
  theme(text = element_text(size=15), 
        axis.text.x = element_text(angle = 90, hjust = 1, face="italic")) 

# ggsave(rank_violin, filename = here("plots", "01_find_tregs", "main_pdf","fig3_rank_violin.png"), width = 3, height = 6)
ggsave(rank_violin, filename = here("plots", "01_find_tregs", "main_pdf","fig3_rank_violin.pdf"), width = 3, height = 6)

rank_violin_ct <- ggplot(rank_long, aes(x = cellType.Broad, y = rank, fill = cellType.Broad)) +
  geom_violin(scale = "width") +
  labs(x = "Cell Type") +
  facet_wrap(~Symbol, ncol = 2) +
  scale_fill_manual(values = cell_colors) +
  theme_bw() +
  theme(text = element_text(size=15), 
        legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(face="italic")) 

# ggsave(rank_violin_ct, filename = here("plots", "01_find_tregs", "main_pdf","fig3_rank_violin_ct.png"), width = 5.5, height = 6)
ggsave(rank_violin_ct, filename = here("plots", "01_find_tregs", "main_pdf","fig3_rank_violin_ct.pdf"), width = 5.5, height = 6)

#### Linear modeling ####
load(here("processed-data", "01_find_tregs", "lmfit.Rdata"), verbose = TRUE) # tt

## bind with rank invar data
invar_t <- tt %>% rownames_to_column("gene") %>%
  left_join(rank_invar_df %>% select(-Symbol)) %>%
  left_join(gene_symbol) %>%
  left_join(genes_of_interest) %>%
  mutate(Signif = cut(adj.P.Val, c(0, 0.001 ,0.01, 0.05, 1), include.lowest = TRUE),
         canidate = Symbol %in% treg_list[1:4],
         dotdot = Symbol %in% dotdotdot_genes,
         probe = Symbol %in% probe_list$Gene.Symbol,
         gene_anno = ifelse(rank_invar %in% 1:5, "RI Worst", gene_anno)
  ) %>%
  replace_na(list(gene_anno = "other"))

invar_t %>% filter(dotdot)
invar_t %>% filter(canidate)
invar_t %>% dplyr::count(gene_anno)
invar_t %>% dplyr::count(probe, gene_anno)
invar_t %>% dplyr::count(probe)
#   probe     n
# 1 FALSE 20238
# 2  TRUE  2800
invar_t %>% dplyr::count(Signif)
#         Signif     n
# 1    [0,0.001] 23000
# 2 (0.001,0.01]    13
# 3  (0.01,0.05]    11
# 4     (0.05,1]    14

invar_t %>% filter(!is.na(rank_invar)) %>% dplyr::count(gene_anno)

invar_t_scatter <- ggplot(invar_t, aes(x = t, y = rank_invar))+
  geom_point(aes(color = probe)) +
  geom_point(data= filter(invar_t, canidate), aes(fill = probe), shape=21, color = "black") +
  geom_text_repel(aes(label = ifelse(canidate, Symbol, NA)), size = 3) +
  labs(title = "Model: log2(counts + 1) ~ log2(sum) + cellType") +
  theme_bw()

ggsave(invar_t_scatter, filename = here("plots", "01_find_tregs", "rank_invar_t_scatter.png"))

pos <- position_jitter(width = 0.3, seed = 2)
invar_t_density <- ggplot(invar_t, aes(x = gene_anno, y = t))+
  geom_jitter(aes(color = Signif), position = pos) +
  geom_text_repel(aes(label = ifelse(gene_anno != "other", Symbol, NA), color = probe), size = 3, position = pos) +
  labs(title = "Model: log2(counts + 1) ~log2(sum) + cellType") +
  theme_bw()

ggsave(invar_t_density , filename = here("plots", "01_find_tregs", "rank_invar_t_densitiy.png"), width = 9)

#### Fig 3. Do Expression trends over cell type track in our HK genes? ####
hk_counts <- log2(as.matrix(assays(sce_pan_unfiltered[genes_of_interest$gene,])$counts)+1)

counts_long <- as.data.frame(hk_counts)%>% 
  rownames_to_column("gene") %>%
  pivot_longer(!gene, names_to = "ID", values_to = "logcount") %>%
  left_join(pd %>% rownames_to_column("ID") %>% select(ID, cellType.Broad,sum)) %>%
  mutate(logsum = log2(sum)) %>%
  left_join(genes_of_interest)%>%
  left_join(invar_t %>% select(gene, t, adj.P.Val))


counts_long$Symbol <- factor(counts_long$Symbol)
counts_long$Symbol <- fct_reorder(counts_long$Symbol, counts_long$t, max)
levels(counts_long$Symbol)

model_anno <- invar_t %>%
  filter(gene_anno != "other") %>%
  mutate(anno = paste("t=", round(t,1),"\nFDR=", scales::scientific(adj.P.Val, didgits = 3)))

model_anno$Symbol <- factor(model_anno$Symbol, levels = levels(counts_long$Symbol))

hk_sum_scatter <- counts_long %>% 
  ggplot(aes(logsum, logcount)) +
  geom_point(aes(color = cellType.Broad), alpha = 0.5, size = 0.5) +
  facet_wrap(~Symbol) +
  geom_text(data = model_anno, aes(label = anno), x = 7.5, y = 12, vjust = "inward", hjust = "inward") +
  scale_color_manual(values = cell_colors) +
  theme_bw() +
  labs(x = "log2(sum)", y = "log2(count + 1)")

ggsave(hk_sum_scatter, filename = here("plots","HK_gene","hk_sum_scatter.png"), width = 12, height = 10)

hk_sum_smooth <- counts_long %>% 
  ggplot() +
  geom_smooth(aes(logsum, logcount, color = cellType.Broad), method = "lm") +
  facet_wrap(~Symbol) +
  geom_text(data = model_anno, aes(label = anno), x = 7.5, y = 12, vjust = "inward", hjust = "inward") +
  scale_color_manual(values = cell_colors) +
  theme_bw() +
  labs(x = "log2(sum)", y = "log2(count + 1)")

ggsave(hk_sum_smooth, filename = here("plots","HK_gene","hk_sum_smooth.png"), width = 12, height = 10)

hk_sum_scatter <- counts_long %>% 
  filter(Symbol == "MALAT1") %>% 
  ggplot(aes(logsum, logcount, color = cellType.Broad)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_smooth(method = "lm") +
  scale_color_manual(values = cell_colors) +
  theme_bw() +
  labs(x = "log2(sum)", y = "log2(count + 1)", title = "MALAT1")

ggsave(hk_sum_scatter, filename = here("plots","HK_gene","hk_sum_scatter_MALAT1.png"), width = 8)

## Figure 3c
hk_sum_scatter_main <- counts_long %>% 
  filter(Symbol == "MALAT1") %>% 
  ggplot(aes(logsum, logcount, color = cellType.Broad)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_smooth(method = "lm") +
  scale_color_manual(values = cell_colors) +
  facet_wrap(~Symbol) +
  theme_bw() +
  coord_equal() +
  labs(x = "log2(sum)", y = "log2(count + 1)") +
  theme(legend.position="None",
        text = element_text(size=15)) 

# ggsave(hk_sum_scatter_main, filename = here("plots","HK_gene","main_pdf","fig3_hk_sum_scatter_MALAT1.png"), width = 4.5, height = 5)
ggsave(hk_sum_scatter_main, filename = here("plots","HK_gene","main_pdf","fig3_hk_sum_scatter_MALAT1.pdf"), width = 4.5, height = 5)


hk_sum_smooth2 <- counts_long %>% 
  filter((gene_anno == "RI Canidate" | Symbol == "POLR2A") &  Symbol != "MALAT1") %>%
  ggplot() +
  geom_smooth(aes(logsum, logcount, color = cellType.Broad), method = "lm") +
  facet_wrap(~Symbol, scales = "free_y", nrow = 1) +
  # geom_text(data = model_anno, aes(label = anno), x = 7.5, y = 12, vjust = "inward", hjust = "inward") +
  scale_color_manual(values = cell_colors) +
  theme_bw() +
  labs(x = "log2(sum)", y = "log2(count + 1)") +
  theme(legend.position = "None")

ggsave(hk_sum_smooth2, filename = here("plots","HK_gene","hk_sum_smooth2.png"), width = 12, height = 5)

## Figure 3d
hk_sum_smooth2_main <- counts_long %>% 
  filter((gene_anno == "RI Canidate" | Symbol == "POLR2A") &  Symbol != "MALAT1") %>%
  ggplot() +
  geom_smooth(aes(logsum, logcount, color = cellType.Broad), method = "lm") +
  facet_wrap(~Symbol, scales = "free_y", ncol = 1) +
  scale_color_manual(values = cell_colors) +
  theme_bw() +
  labs(x = "log2(sum)", y = "log2(count + 1)") +
  theme(legend.position = "None", text = element_text(size=15))

# ggsave(hk_sum_smooth2_main, filename = here("plots","HK_gene","main_pdf","fig3_hk_sum_smooth2.png"), width = 4, height = 5)
ggsave(hk_sum_smooth2_main, filename = here("plots","HK_gene","main_pdf","fig3_hk_sum_smooth2.pdf"), width = 4, height = 5)

gene_slopes <- counts_long %>% group_by(Symbol) %>%
  do(fitQSV = broom::tidy(lm(logcount ~ logsum, data = .))) %>%
  unnest(fitQSV) %>%
  select(Symbol, term, estimate) %>%
  mutate(term = gsub("[^[:alnum:] ]", "", term)) %>%
  pivot_wider(Symbol, names_from = "term", values_from = "estimate") %>%
  column_to_rownames("Symbol")

hk_sum_scatter_AKT3 <- counts_long %>% 
  filter(Symbol == "AKT3") %>%
  ggplot(aes(logsum, logcount, color = cellType.Broad)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_smooth(method = "lm", color = "black") +
  geom_abline(slope = 0.8672058, intercept = -8.7700637, color = "dark red") +
  facet_wrap(~cellType.Broad, scales = "free") +
  scale_color_manual(values = cell_colors) +
  theme_bw()

ggsave(hk_sum_scatter_AKT3, filename = here("plots","HK_gene","hk_sum_scatter_AKT3.png"), width = 11)

pdf(here("plots","HK_gene","hk_sum_scatter_smooth.pdf"), width = 11, height = 8)
for(g in levels(counts_long$Symbol)){
  message(g)
  st = model_anno %>% filter(Symbol == g) %>% pull(anno)
  
  m = gene_slopes[g,'logsum']
  b = gene_slopes[g,'Intercept']
  
  c <- counts_long %>% 
    filter(Symbol == g) %>%
    ggplot(aes(logsum, logcount, color = cellType.Broad))  +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_smooth(method = "lm", color = "black") +
    facet_wrap(~cellType.Broad, scales = "free") +
    geom_abline(slope = m, intercept = b, color = "dark red") +
    scale_color_manual(values = cell_colors) +
    theme_bw() +
    labs(title = g, subtitle = st, x = "log2(sum)", y = "log2(count + 1)")
  
  print(c)
}
dev.off()


# sgejobs::job_single('housekeeper_genes_plots', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript housekeeper_genes_plots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

