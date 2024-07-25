# Create a list of known iLINCS Signatures of the top 3 families

library(creedenzymatic)
library(drugfindR)
library(tidyverse)

targets <- read_csv("results/KRSA_quartile_rankings_aggregated_ordered_all.csv") |>
  slice_head(n = 3) |>
  pull(Kinase)

genes <- kinome_mp_file_v3 |>
  filter(krsa_id %in% targets) |>
  select(krsa_id, hgnc_symbol) |>
  arrange(krsa_id) |>
  write_csv("ancillary/selected_kinase_gene_list.csv")

kd_signatures <- kinome_mp_file_v3 |>
  filter(krsa_id %in% targets) |>
  select(krsa_id, hgnc_symbol) |>
  arrange(krsa_id) |>
  inner_join(drugfindR:::kd_metadata, by = c("hgnc_symbol" = "Source")) |>
  write_csv("ancillary/kd_filter.csv")

oe_signatures <- kinome_mp_file_v3 |>
  filter(krsa_id %in% targets) |>
  select(krsa_id, hgnc_symbol) |>
  arrange(krsa_id) |>
  inner_join(drugfindR:::oe_metadata, by = c("hgnc_symbol" = "Source")) |>
  write_csv("ancillary/oe_filter.csv")
