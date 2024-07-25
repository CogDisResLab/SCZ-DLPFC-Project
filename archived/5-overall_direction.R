# Combined direction of Pairwise Reverse KRSA


library(tidyverse)

dataset <- read_csv("results/filtered_reverse_krsa.csv") |>
  mutate(Direction = case_when(
    LFC >= 0 ~ 1,
    LFC <= -0 ~ -1,
    TRUE ~ 0
  )) |>
  drop_na(LFC)


directions <- dataset |>
  group_by(Kinase, Dataset) |>
  select(Direction) |>
  summarise(Kinase_Direction = sum(Direction)) |>
  summarise(Direction = mean(Kinase_Direction)) |>
  mutate(Direction = round(Direction, 3))

