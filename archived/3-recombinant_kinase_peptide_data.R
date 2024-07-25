# Create Dataset with only the recombinant reporter peptides for each kinase

library(tidyverse)

process_KEA3 <- function(file) {
  col_spec <- cols(Peptide = col_character(),
                   Score = col_double())

  df <- read_delim(file, delim = " ", col_types = col_spec)

  df
}

filter_data <- function(kinase, filtered_peptides, peptide_data) {
  peptides <- filtered_peptides |>
    filter(Source == kinase) |>
    pull(Peptide)

  data <- peptide_data |>
    filter(Peptide %in% peptides)

  data
}

selected_peptides <-
  read_csv("data/selected_peptide_data/top_peptides_per_recombinant_kinase.csv")

peptide_data <- list.files("data/KEA3_Reports/", full.names = TRUE) |>
  set_names(list.files("data/KEA3_Reports/") |> str_extract("R\\dC\\dP\\d")) |>
  map_dfr(process_KEA3, .id = "Pair")

kinases <- selected_peptides |>
  pull(Source) |>
  unique() |>
  set_names()

filtered_peptide_data <- kinases |>
  map(~ filter_data(.x, selected_peptides, peptide_data)) |>
  map2(kinases, ~ write_csv(.x, str_glue("data/selected_peptide_data/{.y}_filtered_peptides_data.csv")))
