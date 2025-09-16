# Creedenzymatic Analysis
# Using shared utility functions

# Source shared functions
source("utils/data_processing_functions.R")

krsa_files <- list.files("results", "krsa", full.names = TRUE)

uka_files <- list.files("results", "uka", full.names = TRUE)

peptide_files <- list.files("results", "dpp", full.names = TRUE)

comparison_names <- krsa_files |>
  basename() |>
  str_extract(".*krsa_table_(.*)\\.csv", 1L)

result <- list(
  krsa_path = krsa_files,
  uka_path = uka_files,
  peptide_path = peptide_files
) |>
  pmap(process_creedenzymatic) |>
  set_names(comparison_names) |>
  imap_dfr(
    ~ write_csv(
      .x,
      file.path("results", str_glue("{.y}_creedenzymatic.csv"))
    ),
    .id = "Comparison"
  )
