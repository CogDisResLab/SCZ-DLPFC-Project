# Generate Activity report for reverse KRSA

library(tidyverse)

dataset <- read_csv("results/filtered_reverse_krsa.csv") |>
  group_by(Kinase) |>
  filter(!is.na(LFC)) |>
  summarise(mean = mean(LFC), sd = sd(LFC))

