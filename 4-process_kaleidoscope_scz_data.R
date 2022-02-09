# Process the SCZ RNASeq data from Kaleidoscope

library(tidyverse)

source_file <- "raw/Kaleidoscope-Data/SCZ_KS_Datasets.rds"

raw_data <- read_rds(source_file)

genes <- read_csv("ancillary/selected_kinase_gene_list.csv") |>
  pull(hgnc_symbol)

clean_data <- raw_data |>
  filter(HGNC_Symbol %in% genes) |>
  select(-ecdfPlot) |>
  write_csv("results/schizophrenia_gene_lookup.csv")
