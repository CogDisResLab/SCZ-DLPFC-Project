# Create kinmap files for Figure 1

library(tidyverse)

pooled_file <- "data/Pooled-Data/UKA_SCZ_P.csv"
region_file <- "data/Region-Data/UKA_SCZ_P.csv"

scale <- 20
color <- "yellow"

pooled_data <- read_csv(pooled_file)
region_data <- read_tsv(region_file)

pooled_clean <- pooled_data |>
  select(`Kinase Name`, `Median Final score`) |>
  rename(Kinase = `Kinase Name`,
         Score = `Median Final score`) |>
  mutate(Size = round(Score) * scale,
         Shape = 0,
         Fill = color,
         Stroke = "black",
         directive = str_glue("@{Shape}:{Size}:{Fill}:{Stroke}:\n{Kinase}"))

region_clean <- region_data |>
  select(`Kinase Name`, `Median Final score`) |>
  rename(Kinase = `Kinase Name`,
         Score = `Median Final score`) |>
  mutate(Size = round(Score) * scale,
         Shape = 0,
         Fill = color,
         Stroke = "black",
         directive = str_glue("@{Shape}:{Size}:{Fill}:{Stroke}\n{Kinase}"))

pooled_kmap <- pooled_clean |>
  pull(directive) |>
  str_c(collapse = "\n\n") |>
  write_file("ancillary/pooled_kinmap_data.txt")

region_kmap <- region_clean |>
  pull(directive) |>
  str_c(collapse = "\n\n") |>
  write_file("ancillary/region_kinmap_data.txt")

