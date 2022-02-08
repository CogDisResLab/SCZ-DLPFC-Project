# Creednzymatic Analysis of Pairwise Data

library(creedenzymatic)
library(tidyverse)

process_KRSA <- function(file) {
  col_spec <- cols(
    .default  = col_double(),
    Kinase = col_character()
  )

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
  col_spec <- cols(
    Peptide = col_character(),
    Score = col_double()
  )

  df <- read_delim(file, delim = " ", col_types = col_spec) %>%
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



krsa_data <- list.files("data/KRSA_Reports/", full.names = TRUE) %>%
  map(process_KRSA) %>%
  set_names(list.files("data/KRSA_Reports/") %>% str_extract("R\\dC\\dP\\d"))

uka_data <- list.files("data/UKA_Reports/", full.names = TRUE) %>%
  map(process_UKA) %>%
  set_names(list.files("data/UKA_Reports/") %>% str_extract("R\\dC\\dP\\d"))

kea3_data <- list.files("data/KEA3_Reports/", full.names = TRUE) %>%
  map(process_KEA3) %>%
  set_names(list.files("data/KEA3_Reports/") %>% str_extract("R\\dC\\dP\\d"))

all_data <- list(krsa_data, uka_data, kea3_data)

common <- all_data %>%
  map(names) %>%
  reduce(intersect)

processed_data <- all_data %>%
  map(~ magrittr::extract(.x, common)) %>%
  set_names(c("KRSA", "UKA", "KEA3")) %>%
  pmap(~ combine_tools(KRSA_df = ..1, UKA_df = ..2, KEA3_df = ..3)) |>
  map2(common, ~ write_csv(.x, file.path("results", "creedenzymatic", str_glue("{.y}_creedenzymatic_results.csv"))))

sig_kinases <- processed_data %>%
  map(~ extract_sig_kinases(.x)) %>%
  map(~ pull(.x, hgnc_symbol)) %>%
  map(~ unique(.x))

final_figures <- processed_data %>%
  map2(sig_kinases, ~ filter(.x, hgnc_symbol %in% {{.y}})) %>%
  map(quartile_figure) %>%
  map2(names(.), ~ ggsave(str_glue("{.y}_Quartile_Figure.png"), plot = .x, path = file.path("figures", "creedenzymatic"), width = 16, height = 8, units = "in"))
