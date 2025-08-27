# Do pairwise correlation of all relevant samples

suppressPackageStartupMessages({
    library(tidyverse)
    library(rstatix)
})

creedencombined_files <- list.files("results", "creedencombined") |>
    set_names(~ str_remove(.x, fixed("_STK_creedencombined.csv")))

combined_data <- creedencombined_files |>
    map(~ read_csv(file.path("results", .x))) |>
    map(~ select(.x, HGNC, Rescaled)) |>
    bind_rows(.id = "Dataset") |>
    pivot_wider(names_from = Dataset, values_from = Rescaled, values_fill = 0L) |>
    select(-HGNC)

correlations <- combined_data |>
    cor_test(method = "spearman") |>
    write_csv(file.path("results", "combined_score_correlations.csv"))
