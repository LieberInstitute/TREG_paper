library("tidyverse")
library("xlsx")
library("here")
library("sessioninfo")

## Supplementary Table 1: Numbers of Different Nuclei over Brain Regions
st1 <- read_csv(here("processed-data","01_find_tregs","supp_tables", "nuclei_counts.csv")) 

cell_types <- c("Astrocyte","Endothelial", "Macrophage", "Microglia", "Mural","Oligodendrocyte",
                "Oligodendrocyte Precursor Cell", "T-cell","Excitatory Neuron", "Inhibitory Neuron")

st1_edit <- st1 %>%
  rename(`Cell Type Abbreviation`  = 1) %>%
  group_by(`Cell Type Abbreviation`) %>%
  mutate(total = sum(AMY, DLPFC, HPC, NAc, sACC)) %>%
  add_column(`Cell Type Name` = cell_types, .before = 1)

## Supplementary Table 2: Detailed gene metrics
st2 <- read_csv(here("processed-data","01_find_tregs","supp_tables", "nuclei_counts.csv")) 

## Supplementary Table 3: RNAscope probe combinations and opal dye assignments

## Supplementary Table 4: Opal dye dilutions used to label probes in RNAscope

## Supplementary Table 5: Scanning protocol exposure times

# Supplementary Table 6: HALO cell counts

