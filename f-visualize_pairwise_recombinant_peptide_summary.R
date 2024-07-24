# Visualize the peptide-level data

library(tidyverse)
library(ggpubr)
library(ggprism)

files <-
  list.files(pattern = "*filtered_peptides_data.csv", recursive = TRUE)

kinases <- str_extract(files, "\\w+_filtered") |>
  str_remove("_filtered") |>
  set_names()

sample_metadata <- read_csv("annotation/sample_matching.csv") |>
  select(Pair = Designation, Gender)

peptide_data <- files |>
  set_names(kinases) |>
  map(~ read_csv(.x, show_col_types = FALSE)) |>
  map(~ inner_join(.x, sample_metadata, by = "Pair"))

peptide_data_summarised <- peptide_data |>
  map(~ summarise(.x, Mean = mean(Score), SD = sd(Score), SEM = sd(Score) / sqrt(n()), .by = c(Gender, Pair)))

kinase_summary_plot <- function(dataset, kinase) {

  pair_order <-
    dataset |> select(Pair, Gender) |> unique() |> arrange(desc(Gender)) |> pull(Pair)

  p <- ggbarplot(dataset, "Pair", "Score", facet.by = "Gender",
   order = pair_order, panel.labs = list(Gender = c("Female", "Male")),
    scales = "free_x", fill = "black", error.plot = "errorbar", add = "mean_se") +
    theme_prism(border = TRUE) +
    theme(axis.text.x = element_text(angle = 45))

  p
}

peptide_data |>
  imap(~ kinase_summary_plot(.x, .y)) |>
  imap(
    ~ ggsave(
      filename = str_glue("{.y}_recombinant_kinase_summary_plot.png"),
      plot = .x,
      path = "figures/recombinant_peptide/summary_plots",
      width = 11,
      height = 8.5,
      units = "in",
      bg = "white"
    )
  )
