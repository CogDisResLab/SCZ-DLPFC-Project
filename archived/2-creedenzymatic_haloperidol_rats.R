# Creedenzymatic for the Rats data

# Process the region level data

library(creedenzymatic)
library(tidyverse)

source_files <- list(
  uka = "data/Rats-Data/Rats_HPD_UKA.csv",

  krsa = "data/Rats-Data/Rats_HPD_KRSA.csv",

  kea3 = "data/Rats-Data/Rats_HPD_KEA3.csv"
)


process_KRSA <- function(file) {
  col_spec <- cols(.default  = col_double(),
                   Kinase = col_character())

  df <- read_csv(file, col_types = col_spec) %>%
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

  df <- read_csv(file, col_types = col_spec) %>%
    select(`Kinase Name`, `Median Final score`) %>%
    rename(Kinase = `Kinase Name`, Score = `Median Final score`) %>%
    read_uka(trns = "abs", sort = "desc")
  df
}

process_KEA3 <- function(file) {
  col_spec <- cols(.default = col_skip(),
                   Peptide = col_character(),
                   LFC = col_double())

  df <- read_csv(file, col_types = col_spec) %>%
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

uka_data <- process_UKA(source_files$uka)

krsa_data <- process_KRSA(source_files$krsa)

kea3_data <- process_KEA3(source_files$kea3)


combined_dataset <- combine_tools(KRSA_df = krsa_data,
                                 UKA_df = uka_data,
                                 KEA3_df = kea3_data) |>
  write_csv(file.path("results", "haloperidol_rats_creedenzymatic_results.csv"))


sig_kinases <- extract_sig_kinases(combined_dataset) |>
  pull(hgnc_symbol)


quartile_fig <- combined_dataset |>
  filter(hgnc_symbol %in% sig_kinases) |>
  quartile_figure() |>
  {\(x) {ggsave(str_glue("HPD_Rats_Quartile_Figure.png"), plot = x, path = file.path("figures", "creedenzymatic"), width = 16, height = 8, units = "in")}}()
