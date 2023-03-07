# Visualize the peptide-level data for specific pairs

library(tidyverse)

files <-
  list.files(pattern = "*filtered_peptides_data.csv", recursive = TRUE)

kinases <- str_extract(files, "\\w+_filtered") |>
  str_remove("_filtered") |>
  set_names()

selected_samples <- c("R1C2P3", "R3C3P5", "R4C1P2", "R1C1P2", "R3C1P1", "R2C2P3")

sample_metadata <- read_csv("annotation/sample_matching.csv") |>
  select(Pair = Designation, Gender) |>
  filter(Pair %in% selected_samples)

peptide_data <- files |>
  set_names(kinases) |>
  map(~ read_csv(.x, show_col_types = FALSE)) |>
  map(~ inner_join(.x, sample_metadata, by = "Pair"))

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

peptide_data |>
  imap(~ comparative_reverse_krsa(.x, .y)) |>
  imap(
    ~ ggsave(
      filename = str_glue("{.y}_subset_recombinant_peptide_activity_boxplot.png"),
      plot = .x,
      path = "figures/recombinant_peptide",
      width = 12,
      height = 8,
      units = "in",
      bg = "white"
    )
  )
