# Similarity between Region and Cell Level datasets

library(tidyverse)

pooled_data <- "data/Pooled-Data/KRSA_SCZ_P.csv" |>
  read_csv() |>
  arrange(Kinase)
region_data <- "data/Region-Data/KRSA_SCZ_P.csv" |>
  read_tsv() |>
  arrange(Kinase)

camk <- read_csv("ancillary/CAMK_gene_subset.csv") |>
  pull(krsa_id) |>
  unique()

camk_pooled <- pooled_data |>
  filter(Kinase %in% camk)

camk_region <- region_data |>
  filter(Kinase %in% camk)

cmgc <- read_csv("ancillary/CMGC_gene_subset.csv") |>
  pull(krsa_id) |>
  unique()

cmgc_pooled <- pooled_data |>
  filter(Kinase %in% cmgc)

cmgc_region <- region_data |>
  filter(Kinase %in% cmgc)

whole_cor <- cor(pooled_data$Z, region_data$Z)

camk_cor <- cor(camk_pooled$Z, camk_region$Z)

cmgc_cor <- cor(cmgc_pooled$Z, cmgc_region$Z)

