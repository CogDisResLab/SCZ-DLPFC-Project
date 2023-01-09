# Coral Input Generator
#

library(tidyverse)

kinome_mp <- creedenzymatic::kinome_mp_file |>
  select(hgnc_symbol, uka_id)

pooled_file <- "data/Pooled-Data/UKA_SCZ_P.csv"
region_file <- "data/Region-Data/UKA_SCZ_P.csv"

pooled_data <- read_csv(pooled_file)
region_data <- read_tsv(region_file)

pooled_clean <- pooled_data |>
  select(`Kinase Name`, `Median Final score`) |>
  rename(Kinase = `Kinase Name`,
         Score = `Median Final score`) |>
  mutate(Kinase = str_to_upper(Kinase)) |>
  left_join(kinome_mp, by = c(Kinase = "uka_id"))

region_clean <- region_data |>
  select(`Kinase Name`, `Median Final score`) |>
  rename(Kinase = `Kinase Name`,
         Score = `Median Final score`) |>
  mutate(Kinase = str_to_upper(Kinase)) |>
  left_join(kinome_mp, by = c(Kinase = "uka_id"))

