# Generate group-wise correlation plot
# Using shared utility functions

# Source shared functions
source("utils/data_processing_functions.R")

labels <- c(
    "C3N_C3NPC", "C3N_D3N", "C3N_D3NPC", "CTL_F_CTL_M", "CTL_SCZ",
    Haloperidol_Control = "Rats", "SCZ_CTL", SCZ_F_CTL_F = "SCZ_F", "SCZ_F_SCZ_M",
    "SCZ_M_CTL_M"
)

# Generate correlation heatmap
p <- generate_correlation_heatmap(
    correlation_file = "results/combined_score_correlations.csv",
    output_file = "groupwise_correlations.png",
    filter_pattern = "SCZ_\\d{1,3}",
    labels = labels,
    title = "Group-wise Correlations"
)
