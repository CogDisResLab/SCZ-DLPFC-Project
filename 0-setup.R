# Setup the directory structure

library(purrr)

create_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}

results_dir <- "results"
data_dir <- "data"
figures_dir <- "figures"
depot_dir <- "depot"


all_dirs <- c(results_dir, data_dir, figures_dir, depot_dir)

all_dirs |>
  walk(create_dir)
