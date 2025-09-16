# Gene Subset Analysis for Kinome Groups
# Consolidated script using shared utility functions

# Source shared functions
source("utils/data_processing_functions.R")

# Define kinome groups to analyze
kinome_groups <- c("CMGC", "CAMK", "AGC", "TK", "TKL", "STE", "CK1", "Other")

# Analyze each group
results <- kinome_groups |>
  set_names() |>
  map(~ analyze_gene_subset(
    group_name = .x,
    kinome_file = "reference_data/kinome_mp_file_v5.csv",
    input_data_file = "data/kaleidoscope_data/KS_SCZ_records.csv",
    output_dir = "results"
  ))

# Optional: Generate summary statistics
summary_stats <- results |>
  map_dfr(~ tibble(
    n_genes = nrow(.x),
    unique_symbols = n_distinct(.x$HGNC_Symbol)
  ), .id = "kinome_group") |>
  write_csv("results/kinome_group_analysis_summary.csv")

message("Gene subset analysis completed for all kinome groups")
print(summary_stats)
