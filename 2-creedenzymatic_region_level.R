# Creedenzymatic analysis on region-level data

# Process the region level data

library(creedenzymatic)
library(tidyverse)

source_files <- list(
  uka = list(
    uka_male = "data/Region-Data/UKA_SCZ_M.csv",
    uka_female = "data/Region-Data/UKA_SCZ_F.csv",
    uka_pooled = "data/Region-Data/UKA_SCZ_P.csv"),

  krsa = list(
    krsa_male = "data/Region-Data/KRSA_SCZ_M.csv",
    krsa_female = "data/Region-Data/KRSA_SCZ_F.csv",
    krsa_pooled = "data/Region-Data/KRSA_SCZ_P.csv"),

  kea3 = list(
    kea3_male = "data/Region-Data/KEA3_SCZ_M.csv",
    kea3_female = "data/Region-Data/KEA3_SCZ_F.csv",
    kea3_pooled = "data/Region-Data/KEA3_SCZ_P.csv")
)


process_KRSA <- function(file) {
  col_spec <- cols(.default  = col_double(),
                   Kinase = col_character())

  df <- read_tsv(file, col_types = col_spec) %>%
    select(Kinase, Z) %>%
    rename(Score = Z) %>%
    read_krsa(trns = "abs", sort = "desc")
  df
}

process_UKA <- function(file) {
  col_spec <- cols(
    .default  = col_double(),
    `Kinase Uniprot ID` = col_character(),
    `Kinase Name` = col_character()
  )

  df <- read_tsv(file, col_types = col_spec) %>%
    select(`Kinase Name`, `Median Final score`) %>%
    rename(Kinase = `Kinase Name`, Score = `Median Final score`) %>%
    read_uka(trns = "abs", sort = "desc")
  df
}

process_KEA3 <- function(file) {
  col_spec <- cols(.default = col_skip(),
                   Peptide = col_character(),
                   LFC = col_double())

  df <- read_tsv(file, col_types = col_spec) %>%
    rename(Score = LFC) %>%
    read_kea(
      filter = T,
      cutoff = 0.5,
      cutoff_abs = T,
      sort = "asc",
      trns = "abs"
    )
  df
}

extract_sig_kinases <- function(df) {
  out <-
    df %>%
    filter(Qrt >= 2) %>%
    select(Kinase, Method, hgnc_symbol) %>%
    mutate(xx = 1) %>%
    pivot_wider(
      names_from = Method,
      values_from = xx,
      values_fill = 0
    ) %>%
    mutate(total = KEA3 + KRSA + UKA) %>%
    filter(total >= 2,
           !is.na(hgnc_symbol))

  out
}

uka_data <- source_files$uka |>
  map(~ process_UKA(.x))

krsa_data <- source_files$krsa |>
  map(~ process_KRSA(.x))


kea3_data <- source_files$kea3 |>
  map(~ process_KEA3(.x))

combined_male <- combine_tools(KRSA_df = krsa_data$krsa_male,
                               UKA_df = uka_data$uka_male,
                               KEA3_df = kea3_data$kea3_male) |>
  write_csv(file.path("results", "region_male_creedenzymatic_results.csv"))

combined_female <- combine_tools(KRSA_df = krsa_data$krsa_female,
                                 UKA_df = uka_data$uka_female,
                                 KEA3_df = kea3_data$kea3_female) |>
  write_csv(file.path("results", "region_female_creedenzymatic_results.csv"))

combined_pooled <- combine_tools(KRSA_df = krsa_data$krsa_pooled,
                                 UKA_df = uka_data$uka_pooled,
                                 KEA3_df = kea3_data$kea3_pooled) |>
  write_csv(file.path("results", "region_pooled_creedenzymatic_results.csv"))

sig_kinases_male <- extract_sig_kinases(combined_male) |>
  pull(hgnc_symbol)
sig_kinases_female <- extract_sig_kinases(combined_female) |>
  pull(hgnc_symbol)
sig_kinases_pooled <- extract_sig_kinases(combined_pooled) |>
  pull(hgnc_symbol)


male_quartile_fig <- combined_male |>
  filter(hgnc_symbol %in% sig_kinases_male) |>
  quartile_figure() |>
  {\(x) {ggsave(str_glue("Male_Region_Quartile_Figure.png"), plot = x, path = file.path("figures", "creedenzymatic"), width = 16, height = 8, units = "in")}}()


female_quartile_fig <- combined_female |>
  filter(hgnc_symbol %in% sig_kinases_female) |>
  quartile_figure() |>
  {\(x) {ggsave(str_glue("Female_Region_Quartile_Figure.png"), plot = x, path = file.path("figures", "creedenzymatic"), width = 16, height = 8, units = "in")}}()

pooled_quartile_fig <- combined_pooled |>
  filter(hgnc_symbol %in% sig_kinases_pooled) |>
  quartile_figure() |>
  {\(x) {ggsave(str_glue("Pooled_Region_Quartile_Figure.png"), plot = x, path = file.path("figures", "creedenzymatic"), width = 16, height = 8, units = "in")}}()
