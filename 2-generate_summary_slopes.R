# Extract the particular peptide values for each pair by diagnosis

library(tidyverse)

files <- list.files("data", "modeled", recursive = TRUE, full.names = TRUE)

pairs_of_interest <- c("R1C2P3", "R3C3P5", "R4C1P2", "R1C1P2", "R3C1P1", "R2C2P3")

sgk1_peptides <- read_csv("data/selected_peptide_data/SGK1_filtered_peptides_data.csv") |>
  pull(Peptide) |>
  unique()

akt1_peptides <- read_csv("data/selected_peptide_data/AKT1_filtered_peptides_data.csv") |>
  pull(Peptide) |>
  unique()

all_data <- files |>
  map_dfr(~read_csv(.x)) |>
  filter(Designation %in% pairs_of_interest)

sgk1_data <- all_data |>
  filter(Peptide %in% sgk1_peptides) |>
  select(SampleName, Group, Peptide, Designation, Gender, slope)

sgk1_summary <- sgk1_data |>
  group_by(Group, Gender) |>
  summarise(AvgSignal = round(mean(slope), 4), SdSignal = round(sd(slope), 4)) |>
  mutate(Result = str_glue("{AvgSignal} ({SdSignal})")) |>
  select(-contains("Signal")) |>
  pivot_wider(names_from = Group, values_from = Result)

sgk1_summary2 <- sgk1_data |>
  group_by(Group) |>
  summarise(AvgSignal = round(mean(slope), 4), SdSignal = round(sd(slope), 4)) |>
  mutate(Result = str_glue("{AvgSignal} ({SdSignal})")) |>
  select(-contains("Signal")) |>
  pivot_wider(names_from = Group, values_from = Result)

akt1_data <- all_data |>
  filter(Peptide %in% akt1_peptides) |>
  select(SampleName, Group, Peptide, Designation, Gender, slope)

akt1_summary <- akt1_data |>
  group_by(Group, Gender) |>
  summarise(AvgSignal = round(mean(slope), 4), SdSignal = round(sd(slope), 4))|>
  mutate(Result = str_glue("{AvgSignal} ({SdSignal})")) |>
  select(-contains("Signal")) |>
  pivot_wider(names_from = Group, values_from = Result)

akt1_summary2 <- akt1_data |>
  group_by(Group) |>
  summarise(AvgSignal = round(mean(slope), 4), SdSignal = round(sd(slope), 4))|>
  mutate(Result = str_glue("{AvgSignal} ({SdSignal})")) |>
  select(-contains("Signal")) |>
  pivot_wider(names_from = Group, values_from = Result)
