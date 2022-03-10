library("tidyverse")
library("here")
library("DeconvoBuddies")
library("SingleCellExperiment")
library("jaffelab")
library("broom")

slope_anno <- function(lm_fit, digits = 1) {
    ci <- confint(lm_fit)
    ## check that coef "cellType.RNAscope.L"
    anno_vals <- c(
        beta = lm_fit$coefficients[[2]],
        ci_low = ci[2, 1],
        ci_high = ci[2, 2]
    )

    anno_vals <- map_chr(anno_vals, ~ as.character(round(.x, digits)))

    slope_anno <- paste0(anno_vals[["beta"]], " (", anno_vals[["ci_low"]], ",", anno_vals[["ci_high"]], ")")
    return(slope_anno)
}

plot_dir <- "plots/02_analyze_halo"
## new cell colors ##
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)
halo_colors <- c(
    cell_colors[c("Oligo", "Excit", "Inhib")],
    c(
        "Multi" = "#4E586A",
        "Other" = "#90A583"
    )
)

halo_colors2 <- halo_colors[c("Oligo", "Excit", "Inhib")]

#### 10x snRNA-seq data ####
load(here("raw-data", "sce_pan.v2.Rdata"), verbose = TRUE)

## filter to just DLPFC data & RI genes
sce_pan <- sce_pan[, sce_pan$region == "dlpfc"]

sce_pan <- sce_pan[rowData(sce_pan)$Symbol %in% c("ARID1B", "AKT3", "MALAT1", "POLR2A"), ]
dim(sce_pan)

sce_pan$cellType.RNAscope <- ifelse(sce_pan$cellType.Broad %in% c("Excit", "Inhib", "Oligo"), as.character(sce_pan$cellType.Broad), "Other")
table(sce_pan$cellType.RNAscope)
table(sce_pan$cellType.Broad)
rownames(sce_pan) <- rowData(sce_pan)$Symbol

nonZero_snRNA <- rowSums(as.matrix(assays(sce_pan)$counts) != 0) / ncol(sce_pan)
# AKT3    ARID1B    MALAT1    POLR2A
# 0.9290305 0.9429566 1.0000000 0.3035172

nonZero_snRNA_cellType <- map_dfc(splitit(sce_pan$cellType.RNAscope), function(i) {
    sce_group <- sce_pan[, i]
    prop_nonzero <- rowSums(as.matrix(assays(sce_group)$counts) != 0) / ncol(sce_group)
    return(prop_nonzero)
})

sn_n_cells <- colData(sce_pan) %>%
    as_tibble() %>%
    group_by(sampleID) %>%
    dplyr::count()

sn_cell_prop <- colData(sce_pan) %>%
    as_tibble() %>%
    group_by(sampleID, cellType.RNAscope) %>%
    summarise(n_cells = n()) %>%
    left_join(sn_n_cells) %>%
    mutate(
        prop = n_cells / n,
        GeneTarget = "snRNA",
        analysis = "snRNA"
    ) %>%
    dplyr::rename(
        Sample = sampleID,
        cell_type = cellType.RNAscope,
        total_cells = n
    ) %>%
    mutate(cell_type = factor(cell_type, levels = names(halo_colors)))

sn_cell_prop$cell_type <- factor(sn_cell_prop$cell_type, levels = names(halo_colors))

#### Sum vs. cell type ####
sum_data <- colData(sce_pan) %>%
    as_tibble() %>%
    select(Sample, uniqueID, cellType.RNAscope, sum) %>%
    mutate(
        cellType.RNAscope = ordered(cellType.RNAscope, levels = c("Excit", "Inhib", "Oligo", "Other")),
        sum_adj = sum / sd(sum)
    ) %>%
    mutate(label = "snRNA-seq")  ## adjust for beta calc

sum_data_main <- sum_data %>%
    filter(cellType.RNAscope != "Other") %>%
    mutate(
        cellType.RNAscope = droplevels(cellType.RNAscope),
        sum_adj = sum / sd(sum)
    )

sd(sum_data$sum)
# [1] 15560.76

sn_fit_main <- lm(sum ~ cellType.RNAscope, data = sum_data_main)
slope_anno(sn_fit_main, 2)
# "-21844.1 (-22172.5,-21515.7)"

## use sum adj
sn_fit_main_adj <- lm(sum_adj ~ cellType.RNAscope, data = sum_data_main)
slope_anno(sn_fit_main_adj, 2)
# "-1.33 (-1.35,-1.31)"

sn_beta_summary <- tibble(RI_gene = "snRNA-seq sum",
                     beta = slope_anno(sn_fit_main, 2),
                     sd = sd(sum_data$sum),
                     beta_adj = slope_anno(sn_fit_main_adj, 2))

## sum over cell type boxplot
sn_sum_boxplot <- sum_data %>%
    ggplot(aes(x = cellType.RNAscope, y = sum, fill = cellType.RNAscope)) +
    geom_boxplot() +
    scale_fill_manual(values = halo_colors) +
    theme_bw() +
    labs(x = "Cell Type") +
    theme(text = element_text(size = 15))

ggsave(sn_sum_boxplot, filename = here(plot_dir, "explore", "sn_sum_boxplot.png"))
## Supp figure 5
ggsave(sn_sum_boxplot + 
           theme(legend.position = "None", axis.text.x = element_text(angle = 90)) +
           facet_wrap(~label),
    filename = here(plot_dir, "supp_pdf", "S5a_sn_sum_boxplot.pdf"), height = 6, width = 3
)

## sum over cell type boxplot
sn_sum_boxplot_main <- sum_data %>%
    filter(cellType.RNAscope != "Other") %>%
    ggplot(aes(x = cellType.RNAscope, y = sum, fill = cellType.RNAscope)) +
    # geom_boxplot() +
    geom_jitter(size = 0.5, alpha = 0.2, color = "grey")+
    geom_boxplot(outlier.shape = NA)+
    scale_fill_manual(values = halo_colors2) +
    facet_wrap(~label) +
    theme_bw() +
    labs(x = "Cell Type") +
    theme(
        text = element_text(size = 15),
        legend.position = "None",
        axis.text.x = element_text(angle = 90)
    )

## Main figure 6
ggsave(sn_sum_boxplot_main,
    filename = here(plot_dir, "main_pdf", "fig6_sn_sum_boxplot.pdf"), height = 6, width = 3
)


#### Load Halo Data ####
load(here("processed-data", "02_analyze_halo", "halo_all.rda"), verbose = TRUE)

halo_all$cell_type <- factor(halo_all$cell_type, names(halo_colors))
## Summary table
cell_type_wide <- halo_all %>%
    filter(RI_gene != "POLR2A", !region_filter) %>%
    group_by(Sample, cell_type) %>%
    dplyr::count() %>%
    pivot_wider(names_from = "cell_type", values_from = "n", names_prefix = "Number ")

(rna_scope_summary <- halo_all %>%
    filter(RI_gene != "POLR2A") %>%
    group_by(Sample) %>%
    summarize(
        `Number cells` = n(),
        `Number cells after filtering` = sum(!region_filter)
    ) %>%
    left_join(cell_type_wide))

write_csv(rna_scope_summary, file = here("processed-data", "02_analyze_halo", "RNAscope_summary.csv"))

## Examine number nuclei
nrow(halo_all)
# [1] 1171800

summary(rna_scope_summary$`Number cells`)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 80796   94528  101475   97691  102654  109130 

summary(rna_scope_summary$`Number cells after filtering`)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 68202   88908   92794   91744   99864  106123

#### Filter Bad Regions ####
halo_all <- halo_all %>% filter(!region_filter)
nrow(halo_all)
# [1] 1099931

#### Number of Cells Per Sample ####
n_cells <- halo_all %>%
    filter(RI_gene != "POLR2A") %>%
    group_by(Sample, GeneTarget) %>%
    summarize(n = n())

n_cells_barplot <- ggplot(n_cells, aes(x = Sample, y = n, fill = GeneTarget)) +
    geom_bar(stat = "identity") +
    labs(y = "Number of Cells") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(n_cells_barplot, filename = here(plot_dir, "explore", "n_cells_sample.png"))

#### Proportion Cells with RI gene ####
RI_hit_barplot <- ggplot(halo_all, aes(x = Sample2, fill = RI_hit)) +
    geom_bar() +
    labs(y = "Number of Cells") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(RI_hit_barplot, filename = here(plot_dir, "explore", "RI_hit_barplot.png"))

n_cells <- n_cells %>% select(Sample, n)

prop_RI_hit <- halo_all %>%
    group_by(RI_gene, Sample) %>%
    filter(RI_hit) %>%
    summarise(n_RI = n()) %>%
    left_join(n_cells) %>%
    mutate(prop_RI = n_RI / n)

mean_n_puncta <- halo_all %>%
    group_by(RI_gene) %>%
    summarize(mean_n_puncta = mean(n_puncta))

(prop_RI_summary <- prop_RI_hit %>%
    group_by(RI_gene) %>%
    summarize(
        mean_prop_RI = mean(prop_RI),
        sd_prop_RI = sd(prop_RI)
    ) %>% left_join(mean_n_puncta))


# A tibble: 4 × 4
# RI_gene mean_prop_RI sd_prop_RI mean_n_puncta
# <chr>          <dbl>      <dbl>         <dbl>
# 1 AKT3           0.881    0.00957          4.09
# 2 ARID1B         0.859    0.00533          3.08
# 3 MALAT1         0.975    0.00293          2.07
# 4 POLR2A         0.785    0.0312           2.75

composition_bar_RI_hit <- prop_RI_summary %>%
    select(RI_gene, mean_prop_TRUE = mean_prop_RI) %>%
    mutate(mean_prop_FALSE = 1 - mean_prop_TRUE) %>%
    pivot_longer(!RI_gene, names_to = "RI_hit", names_prefix = "mean_prop_", values_to = "prop") %>%
    plot_composition_bar(sample_col = "RI_gene", x_col = "RI_gene", ct_col = "RI_hit") +
    theme_bw() +
    labs(fill = "RI Gene Hit")

ggsave(composition_bar_RI_hit, filename = here(plot_dir, "explore", "composition_RI_gene.png"))


#### Cell Type Composition ####
cell_prop <- halo_all %>%
    filter(RI_gene != "POLR2A") %>%
    group_by(Sample, GeneTarget, cell_type, analysis) %>%
    summarise(n_cells = n()) %>%
    left_join(halo_all %>%
        filter(RI_gene != "POLR2A") %>%
        group_by(Sample) %>%
        summarise(total_cells = n())) %>%
    mutate(prop = n_cells / total_cells)

composition_bar_sample <- plot_composition_bar(cell_prop,
    sample_col = "Sample",
    x_col = "Sample",
    min_prop_text = 0.04
) +
    scale_fill_manual(values = halo_colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(composition_bar_sample, filename = here(plot_dir, "explore", "composition_sample.png"), width = 9)

cell_prop2 <- rbind(cell_prop, sn_cell_prop) ## compare with snRNA-seq
cell_prop2$cell_type <- factor(cell_prop2$cell_type, levels = names(halo_colors))

composition_bar_gene <- plot_composition_bar(cell_prop2,
    sample_col = "Sample",
    x_col = "GeneTarget",
    min_prop_text = 0.02
) +
    scale_fill_manual(values = halo_colors) +
    theme_bw()

ggsave(composition_bar_gene, filename = here(plot_dir, "explore", "composition_gene.png"))

## grant version
cell_prop_g <- halo_all %>%
    filter(RI_gene == "AKT3") %>%
    group_by(Sample, cell_type_simple) %>%
    summarise(n_cells = n()) %>%
    left_join(halo_all %>%
        filter(RI_gene != "POLR2A") %>%
        group_by(Sample) %>%
        summarise(total_cells = n())) %>%
    mutate(prop = n_cells / total_cells) %>%
    dplyr::rename(cell_type = cell_type_simple) %>%
    mutate(Sample = gsub("AKT3_Rep#", "Sec ", Sample))

halo_colors_g <- halo_colors[names(halo_colors) != "Multi"]

composition_bar_sample_g <- plot_composition_bar(cell_prop_g,
    sample_col = "Sample",
    x_col = "Sample",
    min_prop_text = 0.04
) +
    scale_fill_manual(values = halo_colors_g) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 17)) +
    labs(y = "Proportion")

ggsave(composition_bar_sample_g, filename = here(plot_dir, "explore", "composition_sample_grant.png"), width = 5, height = 7)

#### Puncta By Cell Type ####
halo_all$cell_type <- ordered(halo_all$cell_type, levels = c("Excit", "Inhib", "Multi", "Oligo", "Other"))

## box plots
n_puncta_boxplot <- halo_all %>%
    ggplot(aes(x = cell_type, fill = cell_type, y = n_puncta)) +
    geom_jitter(size = 0.5, alpha = 0.1, color = "grey")+
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~RI_gene, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = halo_colors) +
    theme_bw() +
    labs(x = "Cell Type", y = "Number Puncta") +
    theme(text = element_text(size = 15), legend.position = "none", axis.text.x = element_text(angle = 90))

ggsave(n_puncta_boxplot, filename = here(plot_dir, "explore", "puncta_boxplot.png"), width = 11)

n_puncta_violin <- halo_all %>%
    ggplot(aes(x = cell_type, fill = cell_type, y = n_puncta)) +
    # geom_jitter(size = 0.5, alpha = 0.1, color = "grey")+
    geom_violin() +
    facet_wrap(~RI_gene, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = halo_colors) +
    theme_bw() +
    labs(x = "Cell Type", y = "Number Puncta") +
    theme(text = element_text(size = 15), legend.position = "none", axis.text.x = element_text(angle = 90))

ggsave(n_puncta_violin, filename = here(plot_dir, "explore", "puncta_violin.png"), width = 11)



## Supp fig 5 version
ggsave(n_puncta_boxplot, filename = here(plot_dir, "supp_pdf", "S5b_puncta_boxplot.pdf"), height = 6)


## Main slopes
# sum ~ cellType.RNAscope
puncta_beta <- halo_all %>%
    filter(cell_type %in% c("Excit", "Inhib", "Oligo")) %>%
    mutate(cell_type = droplevels(cell_type)) %>%
    group_by(RI_gene) %>%
    do(fit = tidy(lm(n_puncta ~ cell_type, data = .), conf.int = TRUE)) %>%
    unnest(fit) %>%
    filter(term == "cell_type.L") %>%
    mutate(beta = paste0(
        round(estimate, 2), " (",
        round(conf.low, 2), ",",
        round(conf.high, 2), ")"
    ))

(beta_summary <- puncta_beta %>%
    select(RI_gene, beta))
# RI_gene beta               
# <chr>   <chr>              
# 1 AKT3    -5.52 (-5.55,-5.49)
# 2 ARID1B  -2.63 (-2.65,-2.6) 
# 3 MALAT1  -1.22 (-1.24,-1.21)
# 4 POLR2A  -3.49 (-3.51,-3.47)

puncta_sd <- halo_all %>%
    group_by(RI_gene) %>%
    filter(cell_type %in% c("Excit", "Inhib", "Oligo")) %>%
    summarise(sd = sd(n_puncta))

## adjust for sd, summarize
beta_summary <- puncta_beta %>%
    left_join(puncta_sd) %>%
    mutate(estimate_adj = estimate/sd,
           conf.low_adj = conf.low/sd,
           conf.high_adj = conf.high/sd,
           beta_adj = paste0(
               round(estimate_adj, 2), " (",
               round(conf.low_adj, 2), ",",
               round(conf.high_adj, 2), ")"
           )) %>%
    select(RI_gene, beta, sd, beta_adj) %>%
    rbind(sn_beta_summary)

# RI_gene       beta                                  sd beta_adj           
# <chr>         <chr>                              <dbl> <chr>              
# 1 AKT3          -5.52 (-5.55,-5.49)                 5.18 -1.07 (-1.07,-1.06)
# 2 ARID1B        -2.63 (-2.65,-2.6)                  3.42 -0.77 (-0.77,-0.76)
# 3 MALAT1        -1.22 (-1.24,-1.21)                 1.53 -0.8 (-0.81,-0.79) 
# 4 POLR2A        -3.49 (-3.51,-3.47)                 3.34 -1.05 (-1.05,-1.04)
# 5 snRNA-seq sum -21844.07 (-22172.45,-21515.68) 15561.   -1.33 (-1.35,-1.31)          

## Main fig 6
n_puncta_boxplot <- halo_all %>%
    filter(cell_type %in% c("Excit", "Inhib", "Oligo"),
           RI_gene != "MALAT1") %>%
    ggplot(aes(x = cell_type, fill = cell_type, y = n_puncta)) +
    geom_jitter(size = 0.5, alpha = 0.1, color = "grey", height = .2)+
    # geom_boxplot()+
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~RI_gene, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = halo_colors2) +
    theme_bw() +
    labs(x = "Cell Type", y = "Number Puncta") +
    theme(text = element_text(size = 15), 
          legend.position = "none", 
          axis.text.x = element_text(angle = 90),
          strip.text.x = element_text(face="italic"))

ggsave(n_puncta_boxplot, filename = here(plot_dir, "main_pdf", "fig6_puncta_boxplot.pdf"), height = 6)

## grant version
n_puncta_boxplot <- halo_all %>%
    filter(!cell_type %in% c("Multi", "Other"), RI_gene == "AKT3") %>%
    ggplot(aes(x = cell_type, fill = cell_type, y = n_puncta)) +
    geom_boxplot() +
    scale_fill_manual(values = halo_colors_g) +
    theme_bw() +
    labs(x = "Cell Type", y = "Number Puncta") +
    theme(text = element_text(size = 20), legend.position = "none")

ggsave(n_puncta_boxplot, filename = here(plot_dir, "explore", "puncta_boxplot_grant.png"))

#### Cell Size by Cell Type ####
nuc_area_boxplot <- halo_all %>%
    filter(RI_gene != "POLR2A") %>% ## Remove duplicate data
    ggplot(aes(x = cell_type, fill = cell_type, y = nucleus_area)) +
    geom_boxplot() +
    scale_fill_manual(values = halo_colors) +
    geom_hline(yintercept = pi * 5^2, color = "dark red") +
    geom_text(label = "Reasonable Area", y = 85, x = "Inhib", color = "dark red") +
    labs(y = "Nucleus Area µm", x = "Cell Type") +
    theme_bw() +
    theme(text = element_text(size = 15), legend.position = "none")

ggsave(nuc_area_boxplot, filename = here(plot_dir, "explore", "nucleus_area_boxplot.png"), width = 9)

## grant version
nuc_area_boxplot <- halo_all %>%
    filter(!cell_type %in% c("Multi", "Other"), RI_gene == "AKT3") %>%
    ggplot(aes(x = cell_type, fill = cell_type, y = nucleus_area)) +
    geom_boxplot() +
    scale_fill_manual(values = halo_colors) +
    labs(y = "Nucleus Area µm", x = "Cell Type") +
    theme_bw() +
    theme(text = element_text(size = 20), legend.position = "none")

ggsave(nuc_area_boxplot, filename = here(plot_dir, "explore", "nucleus_area_boxplot_grant.png"))

## Area vs. puncta ##

puncta_v_size <- halo_all %>%
    ungroup() %>%
    group_by(RI_gene, cell_type) %>%
    do(fitPuncta = tidy(lm(n_puncta ~ nucleus_area, data = .))) %>%
    unnest(fitPuncta) %>%
    filter(term == "nucleus_area") %>%
    mutate(
        p.bonf = p.adjust(p.value, "bonf"),
        p.bonf.sig = p.bonf < 0.05,
        p.bonf.cat = cut(p.bonf,
            breaks = c(1, 0.05, 0.01, 0.005, 0),
            labels = c("<= 0.005", "<= 0.01", "<= 0.05", "> 0.05"),
            include.lowest = TRUE
        ),
        p.fdr = p.adjust(p.value, "fdr"),
        log.p.bonf = -log10(p.bonf)
    )

puncta_v_size_rm0 <- halo_all %>%
    filter(n_puncta != 0) %>%
    group_by(RI_gene, cell_type) %>%
    do(fitPuncta = tidy(lm(n_puncta ~ nucleus_area, data = .))) %>%
    unnest(fitPuncta)

n_puncta_size_scatter <- halo_all %>%
    ggplot(aes(x = nucleus_area, y = n_puncta)) +
    geom_point(aes(color = cell_type), size = 0.2, alpha = 0.2) +
    geom_smooth(method = "lm", color = "black") +
    scale_color_manual(values = halo_colors) +
    facet_grid(RI_gene ~ cell_type) +
    labs(x = "Nucleus Area µm", y = "Number of Puncta") +
    theme_bw() +
    theme(text = element_text(size = 15), legend.position = "none")

ggsave(n_puncta_size_scatter, filename = here(plot_dir, "explore", "puncta_size_scatter.png"), width = 9)

#### Spatial ####
hex_n_cells <- halo_all %>%
    filter(RI_gene != "POLR2A") %>% ## Remove duplicate data & normal nuclei
    ggplot(aes(x = XMax, y = YMax)) +
    geom_hex(bins = 200) +
    scale_fill_continuous(type = "viridis") +
    facet_wrap(~GeneTarget) +
    coord_equal() +
    scale_y_reverse()

ggsave(hex_n_cells, filename = here(plot_dir, "explore", "hex_n_cells.png"))

hex_mean_area <- halo_all %>%
    filter(RI_gene != "POLR2A") %>%
    ggplot(aes(x = XMax, y = YMax, z = nucleus_area)) +
    stat_summary_hex(fun = mean, bins = 100) +
    scale_fill_continuous(type = "viridis", name = "Mean Area") +
    facet_wrap(~GeneTarget, nrow = 1) +
    coord_equal() +
    scale_y_reverse() +
    theme(legend.position = "bottom")

ggsave(hex_mean_area, filename = here(plot_dir, "explore", "hex_mean_area.png"))

#### Check out RI Expression by XY ####
hex_ri <- halo_all %>%
    filter(RI_hit) %>%
    ggplot(aes(x = XMax, y = YMax)) +
    geom_hex(bins = 100) +
    scale_fill_continuous(type = "viridis") +
    facet_wrap(~RI_gene, nrow = 1) +
    coord_equal() +
    scale_y_reverse()

ggsave(hex_ri, filename = here(plot_dir, "explore", "hex_RI.png"), width = 7)

my_breaks <- c(2, 10, 50, 250, 1250, 6000)

hex_ri_sum <- halo_all %>%
    ggplot(aes(x = XMax, y = YMax, z = n_puncta)) +
    stat_summary_hex(fun = sum, bins = 100) +
    scale_fill_continuous(
        type = "viridis", trans = "log", name = "Sum Puncta",
        breaks = my_breaks, labels = my_breaks
    ) +
    facet_wrap(~RI_gene, nrow = 1) +
    coord_equal() +
    scale_y_reverse()

ggsave(hex_ri_sum, filename = here(plot_dir, "explore", "hex_RI_sumlog.png"), width = 9)

hex_ri_sum2 <- halo_all %>%
    ggplot(aes(x = XMax, y = YMax, z = n_puncta)) +
    stat_summary_hex(fun = sum, bins = 100) +
    scale_fill_continuous(type = "viridis", name = "Sum Puncta") +
    facet_wrap(~RI_gene, nrow = 1) +
    coord_equal() +
    scale_y_reverse() +
    theme(legend.position = "bottom")

ggsave(hex_ri_sum2, filename = here(plot_dir, "explore", "hex_RI_sum.png"), width = 9)


hex_ri_mean <- halo_all %>%
    ggplot(aes(x = XMax, y = YMax, z = n_puncta)) +
    stat_summary_hex(fun = mean, bins = 100) +
    scale_fill_continuous(type = "viridis", name = "Mean Puncta") +
    facet_wrap(~RI_gene, nrow = 1) +
    coord_equal() +
    scale_y_reverse() +
    theme(legend.position = "bottom") +
    theme_bw()+
    theme(text = element_text(size = 15), 
          legend.position = "bottom")

ggsave(hex_ri_mean, filename = here(plot_dir, "explore", "hex_mean_puncta.png"), width = 10)
ggsave(hex_ri_mean, filename = here(plot_dir, "supp_pdf", "S5_hex_mean_puncta.pdf"), width = 10)


#### Check out Cell Types by XY ####
walk(levels(halo_all$cell_type), function(ct) {
    ct_hex <- halo_all %>%
        filter(RI_gene != "POLR2A") %>%
        filter(cell_type == ct) %>%
        ggplot(aes(x = XMax, y = YMax)) +
        geom_hex(bins = 100) +
        scale_fill_continuous(type = "viridis") +
        facet_wrap(~GeneTarget, nrow = 1) +
        coord_equal() +
        scale_y_reverse() +
        labs(title = ct)

    ggsave(ct_hex, filename = here(plot_dir, "explore", paste0("hex_ct_", ct, ".png")))
})

#### Figure 4 AKT3 Rep#1 plot ####
halo_a1 <- halo_all %>%
    filter(Sample == "AKT3_Rep#1")

blank_axis <- theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
)


halo_a1_nucleus <- halo_a1 %>%
    ggplot(aes(x = XMax, y = YMax, z = nucleus_area)) +
    stat_summary_hex(fun = mean, bins = 100) +
    scale_fill_continuous(type = "viridis", name = "Mean\nNuclear\nArea") +
    coord_equal() +
    scale_y_reverse() +
    theme_bw() +
    blank_axis

# ggsave(halo_a1_nucleus, filename = here("plots","HK_gene","main_figs","fig5_halo_nucleus.png"),width = 4, height = 2.5)
ggsave(halo_a1_nucleus, filename = here(plot_dir, "main_pdf", "fig5_halo_nucleus.pdf"), width = 4, height = 2.5)

halo_a1_n_puncta <- halo_a1 %>%
    ggplot(aes(x = XMax, y = YMax, z = n_puncta)) +
    stat_summary_hex(fun = mean, bins = 100) +
    scale_fill_continuous(type = "viridis", name = "Mean\nNumber\nPuncta") +
    coord_equal() +
    scale_y_reverse() +
    theme_bw() +
    blank_axis

# ggsave(halo_a1_n_puncta, filename = here("plots","HK_gene","main_figs","fig5_halo_puncta.png"),width = 4, height = 2.5)
ggsave(halo_a1_n_puncta, filename = here(plot_dir, "main_pdf", "fig5_halo_puncta.pdf"), width = 4, height = 2.5)


halo_a1_cell_type <- halo_a1 %>%
    filter(!cell_type %in% c("Multi", "Other")) %>%
    ggplot(aes(x = XMax, y = YMax)) +
    geom_hex(bins = 100) +
    scale_fill_continuous(type = "viridis", name = "Number of Cells") +
    facet_wrap(~cell_type, nrow = 1) +
    coord_equal() +
    scale_y_reverse() +
    theme_bw() +
    blank_axis

# ggsave(halo_a1_cell_type, filename = here("plots","HK_gene","main_figs","fig5_halo_cell_types.png"),width = 8, height = 3)
ggsave(halo_a1_cell_type, filename = here(plot_dir, "main_pdf", "fig5_halo_cell_types.pdf"), width = 8, height = 3)
