# Subset the CAMK dataset from SCZ

library(creedenzymatic)
library(tidyverse)


gene_list <- kinome_mp_file_v5 |>
  filter(group == "CAMK") |>
  select(1:11) |>
  write_csv("ancillary/CAMK_gene_subset.csv") |>
  pull(hgnc_symbol)

filtered_data <- read_csv("data/kaleidoscope_data/KS_SCZ_records.csv") |>
  filter(HGNC_Symbol %in% gene_list) |>
  write_csv("results/CAMK_subset_SCZ_lookup.csv")
