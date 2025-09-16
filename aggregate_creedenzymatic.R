# Reformat and merge creedenzymatic files for interpretation
# Using shared utility functions

# Source shared functions
source("utils/data_processing_functions.R")

creedenzymatic_files <- list.files("results", "creedenzymatic", full.names = TRUE) |>
  set_names(~ str_remove(basename(.x), "_creedenzymatic.csv"))


creedenzymatic_data <- creedenzymatic_files |>
  enframe(value = "file_path") |>
  mutate(
    chip_type = str_extract(name, "[SP]TK"),
    comparison = str_extract(name, "^([A-Z0-9_]+)_[SP]TK", 1L),
    baseline = str_extract(comparison, "([A-Z0-9_]+)_", 1L),
    treatment = str_extract(comparison, "_([A-Z0-9_]+)", 1L),
    data = map(file_path, read_csv),
    pivoted = map(
      data,
      ~ .x |>
        select(Gene = hgnc_symbol, Percentile = Perc, Method) |>
        mutate(
          Method = str_remove(Method, fixed("-")),
          Percentile = round(Percentile * 100L, 2L)
        ) |>
        pivot_wider(
          names_from = Method, values_from = Percentile,
          values_fn = unique, values_fill = -1L
        ) |>
        mutate(Score = combine_score_v(PTMSEA, KEA3, UKA, KRSA)) |>
        select(Gene, Score) |>
        arrange(desc(Score))
    ),
    renamed = map2(
      pivoted, replicate,
      \(x, y) {
        x |>
          rename("replicate_{y}" := Score) # nolint: object_name_linter.
      }
    )
  )

reshaped_creedenzymatic_data <- creedenzymatic_data |>
  select(chip_type, comparison, renamed) |>
  nest(data = -c(chip_type, comparison)) |>
  mutate(
    newdata = map(data, ~ reduce(.x[["renamed"]], full_join, by = "Gene")),
    outfile_name = file.path(
      "results",
      str_glue("{comparison}_{chip_type}_creedenzymatic_aggregated.csv")
    )
  )

walk2(
  reshaped_creedenzymatic_data[["newdata"]],
  reshaped_creedenzymatic_data[["outfile_name"]],
  ~ write_csv(.x, .y)
)
