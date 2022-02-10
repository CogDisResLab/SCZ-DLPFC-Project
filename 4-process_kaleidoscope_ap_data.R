# Process the Kaleidoscope Antipsychotic Data

library(tidyverse)

raw_data <- readRDS("raw/Kaleidoscope-Data/ap_ks_datasets.rds")

genes <- read_csv("ancillary/selected_kinase_gene_list.csv") |>
  pull(hgnc_symbol)

drug_full_name <- function(full_name) {
  drugs = list(
    "CLO" = "Clozapine",
    "HAL" = "Haloperidol",
    "Lox" = "Loxapine",
    "ARI" = "Aripiprazole",
    "QUE" = "Quetiapine",
    "RIS" = "Risperidone",
    "ZIP" = "Ziprasidone",
    "CLZ" = "Clozapine",
    "CL" = "Clozapine",
    "LOX" = "Loxapine",
    "OLA" = "Olanzapine",
    "THI" = "Thiothixene",
    "Valproate" = "Valproate",
    "Lithium" = "Lithium"
  )
  short_name <-
    full_name |>
    str_split(pattern = "[\\W+_]", simplify = TRUE)
  drugs[[short_name[[1]]]]
}

v_drug_name <- Vectorize(drug_full_name)

clean_data <- raw_data |>
  select(-ecdfPlot) |>
  mutate(drug = v_drug_name(DataSet)) |>
  write_csv("data/kaleidoscope_data/KS_AP_records.csv") |>
  filter(HGNC_Symbol %in% genes) |>
  write_csv("results/antipsychotic_gene_lookups.csv")

