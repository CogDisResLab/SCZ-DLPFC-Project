# Summarize the demographic data

library(tidyverse)
library(readxl)

data_file <- "raw/MBC DLPFC Demogs- Kinome Study.xlsx"

sheet <- "Summary"

ctl_data <- read_excel(data_file, sheet, range = "A1:G21", na = c("N/A"))

scz_data <- read_excel(data_file, sheet, range = "I1:O21", na = c("N/A"))

complete_data <- bind_rows(ctl_data, scz_data) |>
  rename(MBC_ID = `MBC#`)

demographics <- complete_data |>
  group_by(DIAG) |>
  summarise(AGEmean = mean(AGE, na.rm = TRUE), AGEsd = sd(AGE, na.rm = TRUE),
            AGEmin = min(AGE, na.rm = TRUE), AGEmax = max(AGE, na.rm = TRUE),
            pHmean = mean(pH, na.rm = TRUE), pHsd = sd(pH, na.rm = TRUE),
            pHmin = min(pH, na.rm = TRUE), pHmax = max(pH, na.rm = TRUE),
            PMImean = mean(PMI, na.rm = TRUE), PMIsd = sd(PMI, na.rm = TRUE),
            PMImin = min(PMI, na.rm = TRUE), PMImax = max(PMI, na.rm = TRUE)) |>
  mutate(across(where(is.numeric), ~ round(.x, 2))) |>
  pivot_longer(cols = -DIAG) |>
  mutate(VAR = str_extract(name, "[pA-Z]+"),
         MEASURE = str_extract(name, "[a-z]{2,}")) |>
  select(-name) |>
  rename(VALUE = value) |>
  pivot_wider(names_from = MEASURE, values_from = VALUE) |>
  mutate(PRINTABLE = str_glue("{mean} +/- {sd} ({min}-{max})")) |>
  write_csv("results/demographics_data_numeric.csv")

sex <- complete_data |>
  group_by(DIAG, SEX) |>
  summarise(COUNT = n()) |>
  pivot_wider(names_from = SEX, values_from = COUNT) |>
  mutate(PRINTABLE = str_glue("{M} Male / {F} Female")) |>
  write_csv("results/demographics_data_sex.csv")



race <- complete_data |>
  group_by(DIAG, RACE) |>
  summarise(COUNT = n()) |>
  pivot_wider(names_from = RACE, values_from = COUNT) |>
  mutate(PRINTABLE = str_glue("{B} Black / {W} White")) |>
  write_csv("results/demographics_data_race.csv")
