library("tidyverse")
library("xlsx")
library("here")

supp_xlxs <- here("processed-data", "03_compile_data", "TREG_SupplementaryTables.xlsx")

## Supplementary Table 1: Numbers of Different Nuclei over Brain Regions
st1 <- read_csv(here("processed-data", "01_find_tregs", "supp_tables", "nuclei_counts.csv"))

cell_types <- c(
    "Astrocyte", "Endothelial", "Macrophage", "Microglia", "Mural", "Oligodendrocyte",
    "Oligodendrocyte Precursor Cell", "T-cell", "Excitatory Neuron", "Inhibitory Neuron"
)

st1_edit <- st1 %>%
    rename(`Cell Type Abbreviation` = 1) %>%
    group_by(`Cell Type Abbreviation`) %>%
    mutate(total = sum(AMY, DLPFC, HPC, NAc, sACC)) %>%
    add_column(`Cell Type Name` = cell_types, .before = 1)

write.xlsx(as.data.frame(st1_edit),
    file = supp_xlxs, sheetName = "Supp Table 1",
    col.names = TRUE, row.names = FALSE, append = FALSE
)

## Supplementary Table 2: Detailed gene metrics
st2 <- read.csv(here("processed-data", "01_find_tregs", "supp_tables", "gene_metrics2.csv"), row.names = 1) %>%
    select(ensembl_id, Symbol, Gene.Type, t)

write.xlsx(as.data.frame(st2),
    file = supp_xlxs, sheetName = "Supp Table 2",
    col.names = TRUE, row.names = FALSE, append = TRUE
)

## Supplementary Table 3: RNAscope probe combinations and opal dye assignments
## Added manually

## Supplementary Table 4: Opal dye dilutions used to label probes in RNAscope
## Added manually

## Supplementary Table 5: Scanning protocol exposure times
## Added manually

# Supplementary Table 6: HALO cell counts
st6 <- read_csv(here("processed-data", "02_analyze_halo", "RNAscope_summary.csv"))

write.xlsx(as.data.frame(st6),
    file = supp_xlxs, sheetName = "Supp Table 6",
    col.names = TRUE, row.names = FALSE, append = TRUE
)
