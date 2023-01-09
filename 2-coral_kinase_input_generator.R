# Create kinmap files for Figure 1

library(tidyverse)

pooled_file <- "data/Pooled-Data/UKA_SCZ_P.csv"
region_file <- "data/Region-Data/UKA_SCZ_P.csv"

pooled_data <- read_csv(pooled_file)
region_data <- read_tsv(region_file)

prepare_coral_data <- function(data, scale = 10) {
  cleaned_data <- data |>
    select(`Kinase Uniprot ID`, `Median Final score`) |>
    rename(Kinase = `Kinase Uniprot ID`,
           Score = `Median Final score`) |>
    mutate(Color = round(Score, 2),
           Size = round(Score, 2))

  size_data <- cleaned_data |>
    select(Kinase, Size)

  color_data <- cleaned_data |>
    select(Kinase, Size)

  out <- list(size = size_data,
              color = color_data)
}

pooled_processed <- prepare_coral_data(pooled_data)
region_processed <- prepare_coral_data(region_data)

