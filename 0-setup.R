# Setup the directory structure

library(purrr)

create_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}

results_dir <- "results"
data_dir <- "data"
krsa_pair_data <- file.path(data_dir, "KRSA_Reports")
uka_pair_data <- file.path(data_dir, "UKA_Reports")
kea3_pair_data <- file.path(data_dir, "KEA3_Reports")
region_data <- file.path(data_dir, "Region-Data")
figures_dir <- "figures"
depot_dir <- "depot"
ancillary_dir <- "ancillary"


all_dirs <- c(results_dir, data_dir, figures_dir, depot_dir,
              krsa_pair_data, kea3_pair_data, uka_pair_data,
              region_data, ancillary_dir)

all_dirs |>
  walk(create_dir)
