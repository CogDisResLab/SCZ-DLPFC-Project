# Significant Peptide Selection by Threshold

library(tidyverse)

files <- list.files("data/KEA3_Reports/", "txt")

paths <- file.path("data", "KEA3_Reports", files)

dataset_names <- str_remove(files, "_KEA3_Peptide_File.txt")

dataset <- paths |>
  set_names(dataset_names) |>
  map(~ read_delim(.x, show_col_types = FALSE, delim = " ")) |>
  imap(~ rename(.x, {
    {
      .y
    }
  } := Score))

all_datasets <- reduce(dataset, full_join, by = "Peptide")

selected_peptides <- all_datasets |>
  pivot_longer(starts_with("R"), names_to = "Pair", values_to = "LogFC") |>
  filter(abs(LogFC) > 0.15) |>
  mutate(value = 1) |>
  select(-LogFC) |>
  pivot_wider(names_from = Pair, values_fill = 0) |>
  mutate(
    Total = R1C1P1 + R1C1P2 + R1C2P3 + R1C3P5 + R2C1P1 + R2C1P2 + R2C2P3 + R2C2P4 + R2C3P5 + R2C3P6 + R3C1P1 + R3C1P2 + R3C2P3 + R3C2P4 + R4C1P1 + R4C1P2 + R1C2P4 + R3C3P5 + R3C3P6,
    across(starts_with("R"), ~ if_else(.x == 0, "No", "Yes"))
  ) |>
  arrange(desc(Total)) |>
  write_csv("results/peptide_selection_by_threshold_0.15.csv")


top_peptides <- all_datasets |>
  pivot_longer(starts_with("R"), names_to = "Pair", values_to = "LogFC") |>
  filter(LogFC > 0.15) |>
  mutate(value = 1) |>
  select(-LogFC) |>
  pivot_wider(names_from = Pair, values_fill = 0) |>
  mutate(
    Total = R1C1P1 + R1C1P2 + R1C2P3 + R1C3P5 + R2C1P1 + R2C1P2 + R2C2P3 + R2C2P4 + R2C3P5 + R2C3P6 + R3C1P1 + R3C1P2 + R3C2P3 + R3C2P4 + R4C1P1 + R4C1P2 + R1C2P4 + R3C3P5 + R3C3P6,
    across(starts_with("R"), ~ if_else(.x == 0, "No", "Yes"))
  ) |>
  arrange(desc(Total)) |>
  write_csv("results/peptide_up_selection_by_threshold_0.15.csv")


bottom_peptides <- all_datasets |>
  pivot_longer(starts_with("R"), names_to = "Pair", values_to = "LogFC") |>
  filter(LogFC < -0.15) |>
  mutate(value = 1) |>
  select(-LogFC) |>
  pivot_wider(names_from = Pair, values_fill = 0) |>
  mutate(
    Total = R2C1P1 + R2C1P2 + R2C3P5 + R3C2P4 + R1C2P4 + R2C2P3 + R3C1P2 + R3C2P3 + R1C1P2 + R1C2P3 + R2C3P6 + R2C2P4 + R3C3P6 + R1C1P1 + R1C3P5 + R3C3P5 + R3C1P1,
    across(starts_with("R"), ~ if_else(.x == 0, "No", "Yes"))
  ) |>
  arrange(desc(Total)) |>
  write_csv("results/peptide_dn_selection_by_threshold_0.15.csv")

all_peptides_signed <- all_datasets |>
  pivot_longer(starts_with("R"), names_to = "Pair", values_to = "LogFC") |>
  filter(abs(LogFC) > 0.15) |>
  mutate(value = case_when(
    LogFC <= -0.15 ~ -1,
    LogFC >= 0.15 ~ 1
  ),
  Direction = if_else(value == -1, "DN", "UP")) |>
  group_by(Peptide, Direction) |>
  mutate(Total = sum(value, na.rm = TRUE)) |>
  pivot_wider(names_from = Direction, values_from = Total) |>
  ungroup() |>
  select(-LogFC) |>
  group_by(Peptide) |>
  fill(UP, DN, .direction = "downup") |>
  pivot_wider(names_from = Pair, values_fill = 0) |>
  mutate(Total_Absolute = sum(abs(UP), abs(DN), na.rm = TRUE),
         Total_Relative = sum(UP, DN, na.rm = TRUE)) |>
  select(Peptide, starts_with("R"), Upregulated = UP, Downregulated = DN, Total_Absolute, Total_Relative) |>
  mutate(across(starts_with("R"), ~ case_when(
    .x == -1 ~ "Down",
    .x == 1 ~ "Up",
    .x == 0 ~ "Absent"
  ))) |>
  write_csv("results/peptide_selection_by_threshold_0.15-enhanced.csv")

