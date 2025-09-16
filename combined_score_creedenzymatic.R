# Calculate the combined score for the creedenzymatic results
# Using shared utility functions

# Source shared functions
source("utils/data_processing_functions.R")

# The main processing logic
creedenzymatic_results <- list.files("results", "creedenzymatic") |>
  set_names(~ str_replace(.x, fixed("creedenzymatic"), "creedencombined"))

results <- creedenzymatic_results |>
  map(process_creeden_file) |>
  imap(
    ~ write_csv(.x, file.path("results", .y))
  ) |>
  imap(~ select(.x, HGNC, CombinedScore = Rescaled)) |>
  bind_rows(.id = "Dataset") |>
  write_csv(file.path("results", "collected_combined_scores.csv"))
