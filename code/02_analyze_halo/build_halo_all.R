library(tidyverse)
library(here)

source("filter_halo.R")
#### load & prep halo data ###
rounds <-paste0("round",1:3)
names(rounds) <- rounds

halo_files <- map(rounds, function(r){
  files <- list.files(path = here("data", "HK_gene", "halo", r), pattern =".csv", full.names = TRUE)
  names(files) <- c("ARID1B", "AKT3", "PM")
  return(files)
})
  
halo_test <- read.csv(halo_files$round2[["ARID1B"]]) ## bad file?

halo_files$round2[["ARID1B"]] <- NA

halo_out <- map_depth(halo_files[c("round1","round3")], 2, read.csv)
halo_colnames <- map_depth(halo_out, 2, colnames)
halo_colnames_unlist <- unlist(halo_colnames)
halo_colnames_edit <- gsub("SLC17A7", "SCL17A7", gsub("POLR2A|MALAT1|ARID1B|AKT3","TREG", gsub("Opal.\\d+","",halo_colnames_unlist)))

hct <- table(unlist(halo_colnames_edit))
length(hct) #75
hct[order(names(hct))]
hct[order(hct)]


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
halo_main %>% dplyr::count(cell_type)

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


halo_all <- rbind(halo_main %>% mutate(analysis = "RNAscope_main"), 
                  halo_PM %>% mutate(GAD1 = NA, analysis = "RNAscope_supp")) %>% 
  separate(Sample, into = c("GeneTarget", "Rep"), sep = "_", remove = FALSE) %>% 
  mutate(Sample2 = paste0(RI_gene, "_", Rep))

halo_all$cell_type <- factor(halo_all$cell_type, levels = c("Oligo", "Excit", "Inhib", "Multi", "Other"))
halo_all$Sample <- factor(halo_all$Sample)

halo_all %>% dplyr::count(Sample)
halo_all %>% dplyr::count(Sample2)
halo_all %>% dplyr::count(cell_type)

#### Print Sample as grid to help filter ####
# load(here("data", "HK_gene", "halo", "halo_all.rda"), verbose = TRUE)
# halo_all <-  halo_all %>% select(-region_filter, -qc_filter)

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
  
  ggsave(grid_hex, filename = here("plots", "HK_gene","halo", paste0("grid_hex_",s,".png")))
})


#### Add filter annotations ####
bad_regions <- read_csv(here("data","HK_gene","halo","bad_regions.csv"))%>%
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
# 1 FALSE         1182371
# 2 TRUE            78040

halo_all <- halo_all %>%
  mutate(problem2 = gsub("[0-9]","",problem)) 

halo_all %>% count(problem2)

#### Plot Regional Filtering ####
anno_test <- ggplot() +
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
  scale_fill_continuous(type = "viridis") +
  facet_wrap(~GeneTarget) +
  scale_x_continuous(minor_breaks = seq(0, 40000, 1000)) +       
  scale_y_continuous(minor_breaks = seq(90000, 0, -1000), trans = "reverse") +
  coord_equal()

ggsave(anno_test, filename = here("plots", "HK_gene","halo", "anno_test.png"), width = 10)


anno_filter_test <- ggplot() +
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
  scale_fill_continuous(type = "viridis") +
  facet_wrap(~GeneTarget) +
  scale_x_continuous(minor_breaks = seq(0, 40000, 1000)) +       
  scale_y_continuous(minor_breaks = seq(90000, 0, -1000), trans = "reverse") +
  coord_equal()

ggsave(anno_filter_test, filename = here("plots", "HK_gene","halo", "anno_filter_test.png"), width = 10)

#### QC box plots ####
## nucleus area
halo_all %>% 
group_by(region_filter) %>%
  summarize(median = median(nucleus_area),
            max = max(nucleus_area))

qc_boxplot_size <- halo_all  %>%
  ggplot(aes(region_filter, nucleus_area, color = problem2)) +
  geom_boxplot()

ggsave(qc_boxplot_size, filename = here("plots", "HK_gene","halo", "qc_boxplot_size.png"))

## n puncta
halo_all %>% group_by(region_filter) %>%
  summarize(median = median(n_puncta),
            max = max(n_puncta))

qc_boxplot_size <- halo_all  %>% 
  ggplot(aes(region_filter, n_puncta, color = problem2)) +
  geom_boxplot()

ggsave(qc_boxplot_size, filename = here("plots", "HK_gene","halo", "qc_boxplot_puncta.png"))

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
save(halo_all, file = here("data", "HK_gene", "halo", "halo_all.rda"))

# sgejobs::job_single('build_halo_all', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript build_halo_all.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
