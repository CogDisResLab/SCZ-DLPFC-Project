# Process the Region-level data

library(tidyverse)

source_files <- list.files(file.path("raw/Region-Data"))

annotation <- source_files |>
  str_match("(\\w{3,4})_(SCZ)_([F|M|P])\\.\\w{3}") |>
  as.data.frame() |>
  set_names(c("file", "app", "disease", "gender")) |>
  as_tibble() |>
  mutate(from = file.path("raw", "Region-Data", file),
         to = file.path("data", "Region-Data", file),
         to = str_replace(to, ".tsv", ".csv"),
         to = str_replace(to, ".txt", ".csv"))

annotation |>
  select(from, to) |>
  pwalk(file.copy)
