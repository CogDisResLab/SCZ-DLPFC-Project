# Calculate the combined score fro the creedenzymatic
# results

suppressPackageStartupMessages({
  library(creedenzymatic)
  library(tidyverse)
})


calculate_weights <- function() {
  kinome <- kinome_mp_file

  total_kinases <- nrow(kinome)
  krsa <- kinome |>
    pull(krsa_id) |>
    keep(~ !is.na(.x)) |>
    length()
  uka <- kinome |>
    pull(uka_id) |>
    keep(~ !is.na(.x)) |>
    length()
  kea3 <- kinome |>
    pull(kea3_id) |>
    keep(~ !is.na(.x)) |>
    length()
  ptmsea <- kinome |>
    pull(ptmsea_id) |>
    keep(~ !is.na(.x)) |>
    length()

  coverage <- c(
    KRSA = round(krsa / total_kinases, 6L),
    UKA = round(uka / total_kinases, 6L),
    KEA3 = round(kea3 / total_kinases, 6L),
    PTMSEA = round(ptmsea / total_kinases, 6L)
  )

  normalized_coverage <- round(coverage / sum(coverage), 4L)

  list(
    coverage = coverage,
    normalized_coverage = normalized_coverage
  )
}

weighted_average_with_coverage <- function(
    p1, p2, p3, p4,
    neutral_value = 0.5,
    penalty_scale = 0.2) {
  calculated_weights_coverage <- calculate_weights()


  calculated_weights <- calculated_weights_coverage$normalized_coverage

  # Compute penalty factors: 1 - (coverage × penalty_scale)
  penalties <- 1L - (calculated_weights_coverage$coverage * penalty_scale)

  # Initialize alpha
  alpha <- 1L

  # Replace missing values (-1) with neutral and apply penalties
  values <- c(p1, p2, p3, p4)
  for (i in 1L:4L) {
    if (values[i] == -1L) {
      alpha <- alpha * penalties[i]
      values[i] <- neutral_value
    }
  }

  # Weighted average
  weighted_avg <- sum(calculated_weights * values)

  # Apply penalty
  final_score <- alpha * weighted_avg
  return(final_score)
}

rescale_combined <- function(x) {
  if (all(x == x[1L])) {
    return(rep(0.5, length(x)))
  } # Handle constant vectors
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


process_creeden_file <- function(filename) {
  filepath <- file.path("results", filename)

  combined_creeden <- read_csv(filepath) |>
    mutate(Method = if_else(Method == "PTM-SEA", "PTMSEA", Method)) |>
    select(Kinase, HGNC = hgnc_symbol, Method, Perc) |>
    mutate(Score = Perc) |>
    pivot_wider(id_cols = HGNC, names_from = Method, values_from = Score, values_fill = -1L, values_fn = unique) |>
    mutate(
      CombinedScore = pmap_dbl(
        list(KRSA, UKA, KEA3, PTMSEA),
        ~ weighted_average_with_coverage(..1, ..2, ..3, ..4)
      ),
      Percentile = ntile(CombinedScore, 100L),
      Rescaled = rescale_combined(CombinedScore)
    )

  combined_creeden
}



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
