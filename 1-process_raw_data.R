# Load the actual raw data and save it in a form useful for plsda

library(tidyverse)
library(KRSA) # nolint: object_name_linter.

read_signal_data <- function(signal_file, sat_file, run) {
  raw_data <- krsa_read(signal_file, sat_file)

  clean_data <- raw_data |>
    dplyr::filter(
      ExposureTime == max(ExposureTime), # nolint: object_usage_linter.
      Cycle == max(Cycle), # nolint: object_usage_linter.
      stringr::str_detect(
        Peptide, "REF", # nolint: object_usage_linter.
        negate = TRUE
      ),
      stringr::str_detect(
        Peptide, "^p",
        negate = TRUE
      ),
      stringr::str_detect(
        SampleName, "^PIM", # nolint: object_usage_linter.
        negate = TRUE
      )
    ) |>
    dplyr::select(Peptide, SampleName, Signal) |> # nolint: object_usage_linter.
    dplyr::mutate(
      Signal = dplyr::if_else(Signal < 0, 0, Signal),
      SampleName = stringr::str_glue("{run}_{SampleName}")
    ) |>
    tidyr::pivot_wider(names_from = SampleName, values_from = Signal)

  clean_data
}

signal_files <- list.files("raw//Paired-Raw-Data",
  pattern = "SigmBg", full.names = TRUE, recursive = TRUE
)
sat_files <- list.files("raw//Paired-Raw-Data",
  pattern = "Sat", full.names = TRUE, recursive = TRUE
)
run <- stringr::str_extract(
  signal_files, "run[0-9]"
)


cleaned_data <- pmap(
  list(signal_files, sat_files, run),
  ~ read_signal_data(..1, ..2, ..3)
) |>
  purrr::reduce(dplyr::full_join, by = "Peptide") |>
  write_csv("peptide_plsda_data.csv")

batch_data <- cleaned_data |>
  pivot_longer(
    cols = -Peptide,
    names_to = "SampleName",
    values_to = "Signal"
  ) |>
  select(SampleName) |>
  unique() |>
  separate_wider_delim(
    cols = SampleName,
    delim = "_",
    names = c("Run", "SampleName"),
    too_many = "merge"
  ) |>
  mutate(SampleName = stringr::str_glue("{Run}_{SampleName}")) |>
  write_csv("batch_plsda_data.csv")
