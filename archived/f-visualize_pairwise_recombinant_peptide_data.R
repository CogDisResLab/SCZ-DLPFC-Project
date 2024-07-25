# Visualize the peptide-level data

library(tidyverse)
library(ggpubr)

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

generate_heatmap <- function(dataset) {
  g <-
    ggplot(data = dataset, aes(x = Pair, y = Peptide, fill = Score))

  p <- g + geom_tile(height = 0.9,
                     width = 0.9,
                     color = "black") +
    theme_minimal() +
    scale_fill_gradient2(
      limits = c(-3.5, 5),
      low = "darkred",
      mid = "white",
      high = "darkgreen"
    ) +
    theme(panel.grid = element_blank())

  p
}

comparative_reverse_krsa <- function(dataset, kinase) {
  g <- ggplot(dataset, aes(x = Pair, y = Score, fill = Gender))

  pair_order <-
    dataset |> select(Pair, Gender) |> unique() |> arrange(desc(Gender)) |> pull(Pair)

  p <- g +
    geom_boxplot(key_glyph = "rect") +
    geom_jitter(width = 0.2,
                height = 0.2,
                show.legend = FALSE) +
    theme_minimal() +
    ggtitle(
      str_glue("Comparison of datasets for {kinase} family"),
      str_glue("High Affinity Recombinant Peptides")
    ) +
    scale_y_continuous(
      name = "Log_2 Fold Change",
      limits = c(-2.0, 2.0),
      breaks = seq(-2.5, 2.5, 0.5)
    ) +
    scale_x_discrete(name = "Dataset", limits = pair_order) +
    scale_fill_manual(
      breaks = c("M", "F"),
      labels = c("Male", "Female"),
      values = c("darkred", "darkgreen")
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) +
    guides(fill = guide_legend())


  p
}

comparative_reverse_krsa_publication <- function(dataset, kinase) {

  pair_order <-
    dataset |> select(Pair, Gender) |> unique() |> arrange(desc(Gender)) |> pull(Pair)

  p <- ggboxplot(dataset, "Pair", "Score", facet.by = "Gender",
   order = pair_order, panel.labs = list(Gender = c("Female", "Male")),
    scales = "free_x", palette = "npg", bxp.errorbar = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}

peptide_data |>
  map(generate_heatmap) |>
  map2(
    kinases,
    ~ ggsave(
      filename = str_glue("{.y}_recombinant_peptide_activity_heatmap.png"),
      plot = .x,
      path = "figures/recombinant_peptide",
      width = 12,
      height = 8,
      units = "in",
      bg = "white"
    )
  )

peptide_data |>
  imap(~ comparative_reverse_krsa(.x, .y)) |>
  imap(
    ~ ggsave(
      filename = str_glue("{.y}_recombinant_peptide_activity_boxplot.png"),
      plot = .x,
      path = "figures/recombinant_peptide",
      width = 12,
      height = 8,
      units = "in",
      bg = "white"
    )
  )

peptide_data |>
  imap(~ comparative_reverse_krsa_publication(.x, .y)) |>
  imap(
    ~ ggsave(
      filename = str_glue("{.y}_recombinant_peptide_activity_boxplot_publication.png"),
      plot = .x,
      path = "figures/recombinant_peptide",
      width = 11,
      height = 8.5,
      units = "in",
      bg = "white"
    )
  )
