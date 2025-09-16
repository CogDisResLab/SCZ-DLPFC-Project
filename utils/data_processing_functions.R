# Data Processing Utility Functions
# Shared functions for kinase analysis across different datasets

# Load required libraries
suppressPackageStartupMessages({
  library(creedenzymatic)
  library(tidyverse)
})

#' Calculate kinome coverage weights
#' @return List containing coverage and normalized coverage weights
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

#' Calculate weighted average with coverage penalty
#' @param p1,p2,p3,p4 Scores for each method
#' @param neutral_value Value to use for missing data
#' @param penalty_scale Scale factor for penalties
#' @return Weighted combined score
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

#' Rescale combined scores to 0-1 range
#' @param x Vector of scores to rescale
#' @return Rescaled vector
rescale_combined <- function(x) {
  if (all(x == x[1L])) {
    return(rep(0.5, length(x)))
  } # Handle constant vectors
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

#' Alternative combine score function with fixed weights and penalties
#' @param ptmsea,kea3,uka,krsa Scores for each method
#' @return Combined score
combine_score <- function(ptmsea, kea3, uka, krsa) {
  # Weights for each method, showing the level of confidence
  # in the method as a percentage.
  weights <- c(0.8, 0.8, 0.9, 1.0)
  # Penalty for each method, showing the level of confidence
  # lost, when the method returns no result.
  penalties <- c(0.9, 0.9, 0.5, 0.6)
  applicable_penalties <- penalties * (c(ptmsea, kea3, uka, krsa) == -1L)
  penalty_factor <- purrr::reduce(
    applicable_penalties[applicable_penalties != 0L],
    `*`,
    .init = 1L
  )

  mean_value <- stats::weighted.mean(c(ptmsea, kea3, uka, krsa),
    weights,
    na.rm = TRUE
  )

  mean_value * penalty_factor
}

#' Vectorized version of combine_score
combine_score_v <- Vectorize(combine_score,
  vectorize.args = c("ptmsea", "kea3", "uka", "krsa")
)

#' Process creedenzymatic analysis with multiple tools
#' @param krsa_path Path to KRSA results file
#' @param uka_path Path to UKA results file
#' @param peptide_path Path to peptide data file
#' @return Combined creedenzymatic results
process_creedenzymatic <- function(krsa_path, uka_path, peptide_path) {
  krsa_data <- readr::read_csv(krsa_path, show_col_types = FALSE) |>
    dplyr::select(Kinase, Score = AvgZ) |>
    creedenzymatic::read_krsa(trns = "abs", sort = "desc")

  uka_data <- readr::read_csv(uka_path, show_col_types = FALSE) |>
    dplyr::select(Kinase = `Kinase Name`, Score = `Median Final score`) |>
    creedenzymatic::read_uka(trns = "abs", sort = "desc")

  peptide_data <- readr::read_csv(peptide_path, show_col_types = FALSE) |>
    dplyr::select(Peptide, Score = totalMeanLFC)

  kea3_data <- creedenzymatic::read_kea(
    peptide_data,
    sort = "asc",
    trns = "abs",
    method = "MeanRank",
    lib = "kinase-substrate"
  )

  ptmsea_data <- creedenzymatic::read_ptmsea(peptide_data)

  combined <- creedenzymatic::combine_tools(
    KRSA_df = krsa_data,
    UKA_df = uka_data,
    KEA3_df = kea3_data,
    PTM_SEA_df = ptmsea_data
  )

  combined
}

#' Perform correlation analysis on combined score data
#' @param results_dir Directory containing creedencombined files
#' @param output_file Output file for correlation results
#' @return Correlation results data frame
perform_correlation_analysis <- function(results_dir = "results", 
                                       output_file = "combined_score_correlations.csv") {
  
  creedencombined_files <- list.files(results_dir, "creedencombined") |>
    set_names(~ stringr::str_remove(.x, stringr::fixed("_STK_creedencombined.csv")))

  combined_data <- creedencombined_files |>
    purrr::map(~ readr::read_csv(file.path(results_dir, .x), show_col_types = FALSE)) |>
    purrr::map(~ dplyr::select(.x, HGNC, Rescaled)) |>
    dplyr::bind_rows(.id = "Dataset") |>
    tidyr::pivot_wider(names_from = Dataset, values_from = Rescaled, values_fill = 0L) |>
    dplyr::select(-HGNC)

  correlations <- combined_data |>
    rstatix::cor_test(method = "spearman") |>
    readr::write_csv(file.path(results_dir, output_file))
  
  correlations
}

#' Generate correlation heatmap
#' @param correlation_file Path to correlation results file
#' @param output_file Output file for the plot
#' @param filter_pattern Pattern to filter out certain comparisons
#' @param labels Named vector of labels for the plot
#' @param title Plot title
#' @param width,height Plot dimensions
#' @return ggplot object
generate_correlation_heatmap <- function(correlation_file, 
                                       output_file,
                                       filter_pattern = "SCZ_\\d{1,3}",
                                       labels = NULL,
                                       title = "Correlation Heatmap",
                                       width = 10L, 
                                       height = 5L) {
  
  correlation_data <- readr::read_csv(correlation_file, show_col_types = FALSE) |>
    dplyr::filter(
      !stringr::str_detect(var1, filter_pattern),
      !stringr::str_detect(var2, filter_pattern)
    ) |>
    dplyr::select(X = var1, Y = var2, Correlation = cor)

  g <- ggplot2::ggplot(correlation_data, ggplot2::aes(x = X, y = Y, fill = Correlation))

  p <- g + 
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::scale_x_discrete(name = NULL, labels = labels, expand = c(0L, 0L)) +
    ggplot2::scale_y_discrete(name = NULL, labels = labels, expand = c(0L, 0L)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(), 
      axis.ticks.y = ggplot2::element_blank(), 
      text = ggplot2::element_text(size = 18L)
    ) +
    ggplot2::ggtitle(title)

  ggplot2::ggsave(output_file,
    width = width,
    height = height,
    units = "in",
    path = "figures", 
    plot = p, 
    bg = "white"
  )
  
  p
}

#' Process KRSA data
#' @param input_data Input data frame containing KRSA results
#' @return Processed data frame with ranks
process_KRSA <- function(input_data) {
  input_data |>
    dplyr::select(Kinase, AvgLog2FC) |>
    dplyr::arrange(desc(AvgLog2FC)) |>
    dplyr::mutate(KRSA_Rank = row_number()) |>
    dplyr::select(Kinase, KRSA_Rank)
}

#' Process UKA data
#' @param input_data Input data frame containing UKA results
#' @return Processed data frame with ranks
process_UKA <- function(input_data) {
  input_data |>
    dplyr::select(Kinase, mLog10P) |>
    dplyr::arrange(desc(mLog10P)) |>
    dplyr::mutate(UKA_Rank = row_number()) |>
    dplyr::select(Kinase, UKA_Rank)
}

#' Process KEA3 data
#' @param input_data Input data frame containing KEA3 results
#' @return Processed data frame with ranks
process_KEA3 <- function(input_data) {
  input_data |>
    dplyr::select(`TF`, `Scaled Rank`) |>
    dplyr::rename(Kinase = `TF`) |>
    dplyr::arrange(`Scaled Rank`) |>
    dplyr::mutate(KEA3_Rank = row_number()) |>
    dplyr::select(Kinase, KEA3_Rank)
}

#' Extract significant kinases from creedenzymatic results
#' @param input_data Input data frame containing creedenzymatic results
#' @param pval_threshold P-value threshold for significance (default: 0.05)
#' @return Vector of significant kinase names
extract_sig_kinases <- function(input_data, pval_threshold = 0.05) {
  input_data |>
    dplyr::filter(log_q_val < -log10(pval_threshold)) |>
    dplyr::pull(kinase_gene_symbol) |>
    unique()
}

#' Process source files with common structure
#' @param source_files List containing file paths for uka, krsa, and kea3
#' @return List of processed data frames
process_source_files <- function(source_files) {
  list(
    uka = readr::read_csv(source_files$uka, show_col_types = FALSE) |> process_UKA(),
    krsa = readr::read_csv(source_files$krsa, show_col_types = FALSE) |> process_KRSA(),
    kea3 = readr::read_csv(source_files$kea3, show_col_types = FALSE) |> process_KEA3()
  )
}

#' Analyze gene subset for kinome groups
#' @param group_name Name of the kinome group (e.g., "CMGC", "CAMK")
#' @param kinome_file Path to kinome mapping file
#' @param input_data_file Path to input data file
#' @param output_dir Directory for output files
#' @return Filtered data frame
analyze_gene_subset <- function(group_name, 
                               kinome_file = "reference_data/kinome_mp_file_v5.csv",
                               input_data_file = "data/kaleidoscope_data/KS_SCZ_records.csv",
                               output_dir = "results") {
  
  # Read kinome mapping file
  kinome_data <- readr::read_csv(kinome_file, show_col_types = FALSE)
  
  # Extract gene list for the specified group
  gene_list <- kinome_data |>
    dplyr::filter(group == group_name) |>
    dplyr::select(1:11) |>
    readr::write_csv(file.path("ancillary", paste0(group_name, "_gene_subset.csv"))) |>
    dplyr::pull(hgnc_symbol)
  
  # Filter input data and save results
  filtered_data <- readr::read_csv(input_data_file, show_col_types = FALSE) |>
    dplyr::filter(HGNC_Symbol %in% gene_list) |>
    readr::write_csv(file.path(output_dir, paste0(group_name, "_subset_SCZ_lookup.csv")))
  
  return(filtered_data)
}

#' Process a single creedenzymatic file
#' @param filename Name of the creedenzymatic file to process
#' @return Processed data frame with combined scores
process_creeden_file <- function(filename) {
  filepath <- file.path("results", filename)

  combined_creeden <- readr::read_csv(filepath, show_col_types = FALSE) |>
    dplyr::mutate(Method = dplyr::if_else(Method == "PTM-SEA", "PTMSEA", Method)) |>
    dplyr::select(Kinase, HGNC = hgnc_symbol, Method, Perc) |>
    dplyr::mutate(Score = Perc) |>
    tidyr::pivot_wider(id_cols = HGNC, names_from = Method, values_from = Score, values_fill = -1L, values_fn = unique) |>
    dplyr::mutate(
      CombinedScore = purrr::pmap_dbl(
        list(KRSA, UKA, KEA3, PTMSEA),
        ~ weighted_average_with_coverage(..1, ..2, ..3, ..4)
      ),
      Percentile = dplyr::ntile(CombinedScore, 100L),
      Rescaled = rescale_combined(CombinedScore)
    )

  combined_creeden
}
