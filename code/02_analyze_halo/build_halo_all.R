library(tidyverse)
library(here)

source("halo_tools.R")
#### load & prep halo data ###
files <- list.files(path = here("raw-data", "HALO", "v3"), pattern =".csv", full.names = TRUE)
names(files) <- c("AKT3", "ARID1B", 'PM')

halo_out <- map(files, read.csv)
map_int(halo_out, nrow)
#   AKT3 ARID1B     PM 
# 301310 285330 292580 

map(halo_out, colnames)

halo_out2 <- map2(halo_out[c("ARID1B", "AKT3")], c("ARID1B", "AKT3"), function(data, gene){
  
  colnames(data) <- gsub(gene, "RI", colnames(data))
  colnames(data) <- gsub("SLC17A7", "SCL17A7", colnames(data))
  
  data2 <- data %>% 
    as_tibble() %>%
    select(Sample = Analysis.Region, 
           ID = Object.Id, 
           MBP = "MBP..Opal.520..Positive", 
           SCL17A7 = "SCL17A7..Opal.690..Positive",
           GAD1 = "GAD1..Opal.620..Positive",
           DAPI.RI,
           cell_area = "Cell.Area..µm..", 
           nucleus_area = "Nucleus.Area..µm..", 
           n_puncta = "RI..Opal.570..Copies",
           XMin, XMax, YMin, YMax) %>%
    mutate(RI_hit = DAPI.RI == 1,
           RI_gene = gene,
           n_marker = MBP + GAD1 + SCL17A7,
           cell_type = case_when(n_marker > 1  ~ "Multi",
                                 GAD1 == 1 ~ "Inhib",
                                 SCL17A7 == 1 ~ "Excit",
                                 MBP == 1 ~ "Oligo",
                                 TRUE~"Other"))
  return(data2)
  
})

halo_main <- do.call("rbind",halo_out2)
halo_main %>% count(cell_type)

halo_PM_copies <- halo_out$PM %>%
  select(Sample = Analysis.Region, 
         ID = Object.Id, 
         n_puncta_MALAT1 = "MALAT1..Opal.520..Copies",
         n_puncta_POLR2A = "POLR2A..Opal.570..Copies") %>%
  pivot_longer(cols = starts_with("n_puncta_"), names_to = "RI_gene",
               names_prefix = "n_puncta_", values_to = "n_puncta")

halo_PM <- halo_out$PM %>% 
  as_tibble() %>%
  select(Sample = Analysis.Region, 
         ID = Object.Id, 
         MBP = "MBP..Opal.620..Positive", 
         SCL17A7 = "SLC17A7..Opal.690..Positive", ##Typo
         DAPI.POLR2A,
         DAPI.MALAT1,
         cell_area = "Cell.Area..µm..", 
         nucleus_area = "Nucleus.Area..µm..",
         XMin, XMax, YMin, YMax) %>%
  pivot_longer(cols = starts_with("DAPI."), names_to = "RI_gene", names_prefix = "DAPI.", values_to = "DAPI.RI") %>%
  mutate(RI_hit = DAPI.RI == 1,
         n_marker = MBP + SCL17A7,
         cell_type = case_when(n_marker > 1 ~ "Multi",
                               SCL17A7 == 1 ~ "Excit",
                               MBP == 1 ~ "Oligo",
                               TRUE~"Other")) %>%
  left_join(halo_PM_copies)

halo_PM %>% dplyr::count(RI_gene, cell_type)


halo_all <- rbind(halo_main %>% mutate(analysis = "RNAscope_1"), 
                  halo_PM %>% mutate(GAD1 = NA, analysis = "RNAscope_2")) %>% 
  separate(Sample, into = c("GeneTarget", "Rep"), sep = "_", remove = FALSE) %>% 
  mutate(Sample2 = paste0(RI_gene, "_", Rep))

halo_all$cell_type <- factor(halo_all$cell_type, levels = c("Oligo", "Excit", "Inhib", "Multi", "Other"))
halo_all$Sample <- factor(halo_all$Sample)

halo_all %>% dplyr::count(Sample)
halo_all %>% dplyr::count(Sample2)
halo_all %>% dplyr::count(cell_type)

#### Print Sample as grid to help filter ####
plot_dir <- "plots/02_analyze_halo"

## print each with grid
walk(levels(halo_all$Sample), function(s){
  message(s)
  halo_s <- halo_all %>% filter(Sample == s, RI_gene != "POLR2A")
  
  grid_hex <- ggplot(halo_s) +
    stat_summary_hex(aes(x = XMax, y = YMax, z = nucleus_area),
                     fun = mean, bins = 100) +
    scale_fill_continuous(type = "viridis") +
    coord_equal() +
    theme_bw() +
    theme(panel.background = element_rect(fill = NA),
          panel.ontop = TRUE,
          panel.grid.major.x = element_line(color="grey60"),
          panel.grid.major.y = element_line(color="grey60")) +  
    scale_x_continuous(minor_breaks = seq(0, 40000, 1000)) +       
    scale_y_continuous(minor_breaks = seq(90000, 0, -1000), trans = "reverse") +
    labs(title = s)
  
  s<- gsub("/","",s)
  
  ggsave(grid_hex, filename = here(plot_dir, "explore", "QC", paste0("grid_hex_",s,".png")))
})


#### Add filter annotations ####
bad_regions <- read_csv(here("processed-data","02_analyze_halo","bad_regions.csv")) %>%
  mutate(GeneTarget = gsub("_Rep#[123]","",Sample))

halo_all <-  halo_all %>% select(-region_filter)
halo_all <- map_dfr(levels(halo_all$Sample), function(s){
  message(s)
  halo_s <- halo_all %>% filter(Sample == s)
  bad_region_s <- bad_regions %>% filter(Sample == s)
  halo_filter <- filter_halo_regions(halo_s, bad_region_s)
  return(halo_filter)
})

halo_all %>% count(region_filter)
# region_filter       n
# 1 FALSE         1099931
# 2 TRUE            71869

halo_all <- halo_all %>%
  mutate(problem2 = gsub("[0-9]","",problem)) 

halo_all %>% count(problem2)

#### Plot Regional Filtering ####
qc_plot <- ggplot() +
  stat_summary_hex(data = halo_all,
                   aes(x = XMax, y = YMax, z = nucleus_area),
                   fun = mean, bins = 100) +
  geom_rect(data = bad_regions,
            aes(xmin = X_min,
                xmax = X_max,
                ymin = Y_min,
                ymax = Y_max),
            fill = NA, color = "red") +
  theme(panel.grid.major.x = element_line(color="grey60"),
        panel.grid.major.y = element_line(color="grey60")) +
  scale_fill_continuous(type = "viridis", name = "Mean\nNuclear\nArea") +
  facet_wrap(~GeneTarget) 

qc_plot_grid <-  qc_plot +
  scale_x_continuous(minor_breaks = seq(0, 40000, 1000)) +       
  scale_y_continuous(minor_breaks = seq(90000, 0, -1000), trans = "reverse") +
  coord_equal()

ggsave(qc_plot_grid, filename = here(plot_dir, "explore", "mean_area_qc_grid.png"), width = 10)

qc_plot_supp <-  qc_plot +
  scale_y_continuous(trans = "reverse") +
  theme_bw() +
  theme(text = element_text(size=15))

ggsave(qc_plot_supp, filename = here(plot_dir, "supp_pdf", "mean_area_qc_grid.png"), width = 10)
ggsave(qc_plot_supp, filename = here(plot_dir, "supp_pdf", "mean_area_qc_grid.pdf"), width = 10)



qc_plot_filtered <- ggplot() +
  stat_summary_hex(data = halo_all %>% filter(!region_filter),
                   aes(x = XMax, y = YMax, z = nucleus_area),
                   fun = mean, bins = 100) +
  geom_rect(data = bad_regions,
            aes(xmin = X_min,
                xmax = X_max,
                ymin = Y_min,
                ymax = Y_max),
            fill = NA, color = "red") +
  scale_fill_continuous(type = "viridis") +
  theme(panel.grid.major.x = element_line(color="grey60"),
        panel.grid.major.y = element_line(color="grey60")) +
  scale_fill_continuous(type = "viridis", name = "Mean\nNuclear\nArea") +
  facet_wrap(~GeneTarget) +
  scale_x_continuous(minor_breaks = seq(0, 40000, 1000)) +       
  scale_y_continuous(minor_breaks = seq(90000, 0, -1000), trans = "reverse") +
  coord_equal()

ggsave(qc_plot_filtered, filename = here(plot_dir, "explore", "mean_area_qc_grid_filter.png"), width = 10)

#### QC box plots ####
## nucleus area
halo_all %>% 
group_by(region_filter) %>%
  summarize(median = median(nucleus_area),
            max = max(nucleus_area))

qc_boxplot_size <- halo_all  %>%
  ggplot(aes(region_filter, nucleus_area, color = problem2)) +
  geom_boxplot()

ggsave(qc_boxplot_size, filename = here(plot_dir, "explore", "qc_boxplot_size.png"))

## n puncta
halo_all %>% group_by(region_filter) %>%
  summarize(median = median(n_puncta),
            max = max(n_puncta))

qc_boxplot_size <- halo_all  %>% 
  ggplot(aes(region_filter, n_puncta, color = problem2)) +
  geom_boxplot()

ggsave(qc_boxplot_size, filename = here(plot_dir, "explore", "qc_boxplot_puncta.png"))

## filter summary
halo_all %>%
  filter(RI_gene != "POLR2A") %>%
  group_by(Sample) %>%
  summarize(n = n(),
            filtered = sum(region_filter),
            prop = filtered/n)

## what cell types end up filtered?
halo_all  %>%
  filter(RI_gene != "POLR2A") %>%
  group_by(cell_type) %>%
  summarise(n_cells = n(),
            n_cells_filtered = sum(region_filter)) %>%
  mutate(prop_filtered = n_cells_filtered/n_cells)


#### Save Data ####
save(halo_all, file = here("processed-data","02_analyze_halo", "halo_all.rda"))

# sgejobs::job_single('build_halo_all', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript build_halo_all.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
