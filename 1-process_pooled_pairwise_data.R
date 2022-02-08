# Process the Pooled Pairwise data

library(tidyverse)

source_path <- file.path("raw", "pooled-Pairwise_data")
source_files <- list.files(source_path, "txt")
dest_path <- file.path("data", "Pooled-Data")

annotation <- source_files |>
  as_tibble() |>
  rename(old_file = value) |>
  mutate(new_file = case_when(
    str_detect(old_file, "KEA3") ~ "KEA3_SCZ_P.csv",
    str_detect(old_file, "KRSA") ~ "KRSA_SCZ_P.csv",
    str_detect(old_file, "Summary") ~ "UKA_SCZ_P.csv"
  ),
  data = map(file.path(source_path, old_file), ~ read_tsv(.x)),
  out = map2(data, file.path(dest_path, new_file), ~ write_csv(.x, .y)))
