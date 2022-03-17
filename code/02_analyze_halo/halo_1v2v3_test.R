library(tidyverse)
library(jaffelab)
library(here)

source("halo_tools.R")
#### load & prep halo data ###
files <- list.files(path = here("raw-data", "HALO", "v3"), pattern = ".csv", full.names = TRUE)
files <- list.files(path = here("raw-data", "HALO", "v1"), pattern = ".csv", full.names = TRUE)

v1 <- read.csv("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/TREG_paper/raw-data/HALO/v1/11112021_RNAscopeInvariance_Arid1b_Fused.tif_object_Data.csv")
v2 <- read.csv("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/TREG_paper/raw-data/HALO/v2/2022-02-18_TREG_ARID1B_Rerun_object_Data_Cyt_Nuc.csv")
v3 <- read.csv("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/TREG_paper/raw-data/HALO/v3/2022-02-18_TREG_ARID1B_Rerun_object_Data_Nuc.csv")

## fix colnames
colnames(v1) <- gsub('_$', '', gsub('(\\.)+', '_', colnames(v1)))
colnames(v2) <- gsub('_$', '', gsub('(\\.)+', '_', colnames(v2)))
colnames(v3) <- gsub('_$', '', gsub('(\\.)+', '_', colnames(v3)))

# Cytoplasm cols are unique to 
colnames(v2)[!colnames(v2) %in% colnames(v3)]
# [1] "MBP_Opal_520_Positive_Cytoplasm"      "MBP_Opal_520_Cytoplasm_Intensity"     "SLC17A7_Opal_690_Positive_Cytoplasm" 
# [4] "SLC17A7_Opal_690_Cytoplasm_Intensity" "GAD1_Opal_620_Positive_Cytoplasm"     "GAD1_Opal_620_Cytoplasm_Intensity"   
# [7] "DAPI_DAPI_Positive_Cytoplasm"         "DAPI_DAPI_Cytoplasm_Intensity"        "Cytoplasm_Area_µm"  

# All V2 colanmes are in V3
colnames(v3)[!colnames(v3) %in% colnames(v2)]
# character(0)

## same n rows between v1 & 2
dim(v1)
# [1] 304505     42
dim(v2)
# [1] 304505     45

## different number of rows in v3
dim(v3)
# [1] 285330     36

v1 <- as_tibble(v1)
v2 <- as_tibble(v2)
v3 <- as_tibble(v3)

## can we merge these?
v1 %>% count(Analysis_Region)
# Analysis_Region      n
# <chr>            <int>
# 1 ARID1B_Rep#1    107334
# 2 ARID1B_Rep#2    101761
# 3 ARID1B_Rep#3     95410

v2 %>% count(Analysis_Region)
# Analysis_Region      n
# <chr>            <int>
# 1 ARID1B_Rep#1    107334
# 2 ARID1B_Rep#2    101761
# 3 ARID1B_Rep#3     95410
v3 %>% count(Analysis_Region)
# Analysis_Region      n
# <chr>            <int>
# 1 ARID1B_Rep#1    101513
# 2 ARID1B_Rep#2     94528
# 3 ARID1B_Rep#3     89289

test_object12 <- v1 %>%
  select(2:8) %>%
  # head(100) %>% 
  rename_with(.fn = ~ paste0(.x, "_v2"), .cols = !c(Object_Id, Analysis_Region)) %>%
  left_join(v2 %>% select(2:8)) %>%
  mutate(d_min = XMin_v2 - XMin,
         d_max = XMax_v2 - XMax,
         q_min = XMin_v2/XMin) 

 test_object23 <- v2 %>%
  select(2:8) %>%
  # head(100) %>% 
  rename_with(.fn = ~ paste0(.x, "_v2"), .cols = !c(Object_Id, Analysis_Region)) %>%
  left_join(v3 %>% select(2:8)) %>%
   mutate(d_min = XMin_v2 - XMin,
          d_max = XMax_v2 - XMax,
          q_min = XMin_v2/XMin) 
 
# Not consistantly different 
test_object  %>%
   count(d_min) %>% 
   arrange(-n)
 
## test joining by X Y 
test_xy <- v2 %>%
  select(2:8) %>%
  # head(100) %>%
  rename_with(.fn = ~ paste0(.x, "_v2"), .cols = c(Object_Id, Algorithm_Name)) %>%
  left_join(v3 %>% select(2:8))

## no coords match :(
test_xy %>%
  count(is.na(Object_Id))
# `is.na(Object_Id)`      n
# <lgl>               <int>
#   1 TRUE               304505


## how do nuc vs nuc + cyto change # puncta
halo_v2 <- v2 %>%
  select(
    Sample = Analysis_Region,
    ID = Object_Id,
    MBP = "MBP_Opal_520_Positive",
    SCL17A7 = "SLC17A7_Opal_690_Positive",
    GAD1 = "GAD1_Opal_620_Positive",
    DAPI_RI = "DAPI_ARID1B",
    cell_area = "Cell_Area_µm",
    nucleus_area = "Nucleus_Area_µm",
    cytoplasm_area = "Cytoplasm_Area_µm",
    cell_area = "Cell_Area_µm",
    n_puncta = "ARID1B_Opal_570_Copies",
    XMin, XMax, YMin, YMax
  ) %>%
  mutate(
    RI_hit = DAPI_RI == 1,
    RI_gene = "ARID1B",
    n_marker = MBP + GAD1 + SCL17A7,
    cell_type = case_when(
      n_marker > 1 ~ "Multi",
      GAD1 == 1 ~ "Inhib",
      SCL17A7 == 1 ~ "Excit",
      MBP == 1 ~ "Oligo",
      TRUE ~ "Other"
    ))

halo_v2 %>%
  count(cell_type)
# cell_type      n
# <chr>      <int>
# 1 Excit      68426
# 2 Inhib      15475
# 3 Multi      11033
# 4 Oligo      51169
# 5 Other     158402

## test plot
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

n_puncta_boxplot <- halo_v2 %>%
  filter(cell_type %in% c("Excit", "Inhib", "Oligo")) %>%
  ggplot(aes(x = cell_type, fill = cell_type, y = n_puncta)) +
  geom_boxplot()+
  # geom_jitter(size = 0.5, alpha = 0.1, color = "grey", height = .2)+
  # geom_boxplot(outlier.shape = NA) +
  facet_wrap(~RI_gene, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = halo_colors2) +
  theme_bw() +
  labs(x = "Cell Type", y = "Number Puncta", title = "HALO v2", subtitle = "nuc + cyto") +
  theme(text = element_text(size = 15), 
        legend.position = "none", 
        axis.text.x = element_text(angle = 90),
        strip.text.x = element_text(face="italic"))

ggsave(n_puncta_boxplot, filename = here(plot_dir, "explore", "puncta_boxplot_ARID1B_v2.png"), height = 6, width = 3)

area_scatter <- halo_v2 %>%
  ggplot(aes(x = nucleus_area, y = cell_area, color = cell_type)) +
  geom_point(alpha = 0.2) + 
  facet_wrap(~cell_type) +
  scale_color_manual(values = halo_colors) +
  theme_bw()

ggsave(area_scatter, filename = here(plot_dir, "explore", "area_scatter_ARID1B_v2.png"))


