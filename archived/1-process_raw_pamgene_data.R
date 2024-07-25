# Process the raw files to generate slopes.
#

library(tidyverse)
library(KRSA)

metadata <- read_csv("annotation/sample_matching.csv") |>
  select(Designation, CTL, SCZ, Gender) |>
  pivot_longer(cols = c(SCZ, CTL), values_to = "SampleName") |>
  select(-name)

signal_files <- list.files("raw", "SigmBg", full.names = TRUE, recursive = TRUE) |>
  keep(~ str_detect(.x, "Paired"))
saturation_files <- list.files("raw", "Saturation", full.names = TRUE, recursive = TRUE) |>
  keep(~ str_detect(.x, "Paired"))

processed_data <- map2(signal_files, saturation_files, krsa_read) |>
  map(~ filter(.x, str_detect(SampleName, "SCZ") | str_detect(SampleName, "CTL"))) |>
  map(krsa_qc_steps) |>
  map(~ mutate(.x, Group = str_extract(SampleName, "^\\w{3}"))) |>
  set_names(str_c("run", 1:4))

processed_end_maxexp <- processed_data |>
  map(~ krsa_extractEndPointMaxExp(.x, "STK"))

processed_end <- processed_data |>
  map(~ krsa_extractEndPoint(.x, "STK"))

pass_peps <- processed_end_maxexp |>
  map(~ krsa_filter_lowPeps(.x, 5))

modeled <- processed_end |>
  map2(pass_peps, ~ krsa_scaleModel(.x, .y)) |>
  map(~ pluck(.x, "scaled")) |>
  map(~ inner_join(.x, metadata, by = "SampleName")) |>
  imap(~ write_csv(.x, str_glue("data/Pairwise-Data/{.y}-modeled-scaled.csv")))



