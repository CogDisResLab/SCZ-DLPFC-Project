# Consolidated Correlation Analysis and Plotting
# Using shared utility functions

# Source shared functions
source("utils/data_processing_functions.R")

# Perform correlation analysis
message("Performing correlation analysis...")
correlations <- perform_correlation_analysis()

# Define labels for plotting
labels <- c(
  "C3N_C3NPC", "C3N_D3N", "C3N_D3NPC", "CTL_F_CTL_M", "CTL_SCZ",
  Haloperidol_Control = "Rats", "SCZ_CTL", SCZ_F_CTL_F = "SCZ_F", 
  "SCZ_F_SCZ_M", "SCZ_M_CTL_M"
)

# Generate group-wise correlation heatmap
message("Generating group-wise correlation heatmap...")
groupwise_plot <- generate_correlation_heatmap(
  correlation_file = "results/combined_score_correlations.csv",
  output_file = "groupwise_correlations.png",
  filter_pattern = "SCZ_\\d{1,3}",
  labels = labels,
  title = "Group-wise Correlations",
  width = 10L,
  height = 5L
)

# Generate pairwise correlation heatmap (including SCZ samples)
message("Generating pairwise correlation heatmap...")
pairwise_plot <- generate_correlation_heatmap(
  correlation_file = "results/combined_score_correlations.csv",
  output_file = "pairwise_correlations.png",
  filter_pattern = NULL, # Include all comparisons
  labels = NULL, # Use original names
  title = "Pairwise Correlations",
  width = 12L,
  height = 8L
)

message("Correlation analysis and plotting completed!")
message("Output files:")
message("- results/combined_score_correlations.csv")
message("- figures/groupwise_correlations.png") 
message("- figures/pairwise_correlations.png")
