# Pathway Analysis - enrichr

library(tidyverse)
library(enrichR)
library(writexl)

stk_id <- readRDS("reference_data/stk_id_map.Rds")
stk_hgnc <- readRDS("reference_data/stk_hgnc_map.Rds")

stk_map <- stk_id |>
  inner_join(stk_hgnc)

data <- read_csv("results/peptide_selection_by_threshold_0.15-enhanced.csv") |>
  select(-starts_with("R")) |>
  mutate(Rank = rank(desc(Total_Relative), ties.method = "max"))

genes_up <- data |>
  filter(Rank <= 6) |>
  select(Peptide) |>
  inner_join(stk_map) |>
  pull(Gene)

genes_dn <- data |>
  filter(Rank >= 92) |>
  select(Peptide) |>
  inner_join(stk_map) |>
  pull(Gene)

genes_all <- c(genes_up, genes_dn)

dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021",
         "GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2022")

up_enriched <- enrichr(genes_up, dbs) |>
  imap(~ write_csv(.x, str_glue("results/Top-{.y}-Pathways.csv")))

writexl::write_xlsx(up_enriched, "results/Up-Pathways.xlsx")

dn_enriched <- enrichr(genes_up, dbs) |>
  imap(~ write_csv(.x, str_glue("results/Bottom-{.y}-Pathways.csv")))

writexl::write_xlsx(up_enriched, "results/Down-Pathways.xlsx")

all_enriched <- enrichr(genes_up, dbs) |>
  imap(~ write_csv(.x, str_glue("results/All-{.y}-Pathways.csv")))

writexl::write_xlsx(up_enriched, "results/All-Pathways.xlsx")
