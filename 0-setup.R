# Setup the directory structure

library(purrr)

create_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}

results_dir <- "results"
creedenzymatic_results_dir <- file.path(results_dir, "creedenzymatic")
data_dir <- "data"
krsa_pair_data <- file.path(data_dir, "KRSA_Reports")
uka_pair_data <- file.path(data_dir, "UKA_Reports")
kea3_pair_data <- file.path(data_dir, "KEA3_Reports")
region_data <- file.path(data_dir, "Region-Data")
rats_data <- file.path(data_dir, "Rats-Data")
pooled_pairwise_data <- file.path(data_dir, "Pooled-Data")
kaleidoscope_data_dir <- file.path(data_dir, "kaleidoscope_data")
figures_dir <- "figures"
depot_dir <- "depot"
ancillary_dir <- "ancillary"
creedenzymatic_figure_dir <- file.path(figures_dir, "creedenzymatic")
reverse_krsa_figure_dir <- file.path(figures_dir, "reverse_krsa")


all_dirs <- c(results_dir, data_dir, figures_dir, depot_dir,
              krsa_pair_data, kea3_pair_data, uka_pair_data,
              region_data, ancillary_dir, creedenzymatic_figure_dir, creedenzymatic_results_dir, pooled_pairwise_data, reverse_krsa_figure_dir,
              kaleidoscope_data_dir, rats_data)

all_dirs |>
  walk(create_dir)
