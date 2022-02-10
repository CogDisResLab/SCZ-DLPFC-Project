# Analyze the CAMK vs CMGC Significance dataset

library(tidyverse)

camk <- read_csv("results/CAMK_subset_SCZ_lookup.csv") |>
  mutate(indicator = case_when(
    P_Value < 0.05 & Log2FC > 1 ~ "HiSig",
    P_Value < 0.05 ~ "Sig",
    Log2FC > 1 ~ "Hi",
    TRUE ~ "None"
  ),
  Family = "CAMK")

cmgc <- read_csv("results/CMGC_subset_SCZ_lookup.csv") |>
  mutate(indicator = case_when(
    P_Value < 0.05 & Log2FC > 1 ~ "HiSig",
    P_Value < 0.05 ~ "Sig",
    Log2FC > 1 ~ "Hi",
    TRUE ~ "None"
  ),
  Family = "CMGC")

combined <- bind_rows(camk, cmgc) |>
  group_by(HGNC_Symbol, indicator, Family) |>
  count() |>
  pivot_wider(names_from = indicator, values_from = n, values_fill = 0) |>
  write_csv("results/Kinase_Familywise_Transcriptional_Activity.csv")
