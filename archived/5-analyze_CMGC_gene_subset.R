# Subset the CMGC dataset from SCZ

library(creedenzymatic)
library(tidyverse)


gene_list <- kinome_mp_file_v5 |>
  filter(group == "CMGC") |>
  select(1:11) |>
  write_csv("ancillary/CMGC_gene_subset.csv") |>
  pull(hgnc_symbol)

filtered_data <- read_csv("data/kaleidoscope_data/KS_SCZ_records.csv") |>
  filter(HGNC_Symbol %in% gene_list) |>
  write_csv("results/CMGC_subset_SCZ_lookup.csv")
