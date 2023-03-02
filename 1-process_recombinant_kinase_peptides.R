# Process reporter peptides for our hit kinases

library(tidyverse)

top_pep_files <-
  list.files(pattern = "*_top_pep[ts]\\w*.csv", recursive = TRUE)

filter_top_hits <- function(df) {
  if ("Cluster" %in% names(df)) {
    df |>
      filter(Cluster == "High Affinity")
  } else {
    df
  }
}

process_top_peptides <- function(file) {
  gene_name <- str_extract(file, "(\\w*)_top") |>
    str_remove("_top")
  dataset <- read_csv(file) |>
    filter_top_hits() |>
    select(-any_of(c("...1", "Threshold", "LFC Threshold", "Cluster", "Type", "ID", "Gene"))) |>
    mutate(Source = gene_name)

  dataset
}

final_top_peptides <- top_pep_files |>
  map_dfr(process_top_peptides) |>
  select(-X) |>
  mutate(Source = if_else(str_detect(Source, "ERK5"), "ERK5", Source)) |>
  unique() |>
  write_csv("data/selected_peptide_data/top_peptides_per_recombinant_kinase.csv")
