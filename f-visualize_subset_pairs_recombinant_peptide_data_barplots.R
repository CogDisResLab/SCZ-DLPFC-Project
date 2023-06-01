# Visualize the peptide-level data for specific pairs as barplots

library(tidyverse)

files <-
  list.files(pattern = "*filtered_peptides_data.csv", recursive = TRUE)

kinases <- str_extract(files, "\\w+_filtered") |>
  str_remove("_filtered") |>
  set_names()

selected_samples <- c("R1C2P3", "R2C2P4", "R3C3P5", "R4C1P2", "R1C1P2", "R2C2P3", "R3C1P1", "R3C2P3")

sample_metadata <- read_csv("annotation/sample_matching.csv") |>
  select(Pair = Designation, Gender) |>
  filter(Pair %in% selected_samples)

peptide_data <- files |>
  set_names(kinases) |>
  map(~ read_csv(.x, show_col_types = FALSE)) |>
  map(~ inner_join(.x, sample_metadata, by = "Pair"))

comparative_reverse_krsa <- function(dataset, kinase) {
  g <- ggplot(dataset, aes(x = Pair, y = Score, fill = Gender))

  p <- g +
    geom_bar(key_glyph = "rect", stat = "summary", fun = "mean") +
    theme_minimal() +
    ggtitle(
      str_glue("Comparison of datasets for {kinase} family"),
      str_glue("High Affinity Recombinant Peptides")
    ) +
    scale_y_continuous(
      name = "Log_2 Fold Change",
      limits = c(-0.25, 0.75),
      breaks = seq(-2.5, 2.5, 0.25)
    ) +
    scale_fill_manual(
      breaks = c("M", "F"),
      labels = c("Male", "Female"),
      values = c("darkred", "darkgreen")
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.direction = "horizontal",
      strip.background = element_rect(fill = "grey80"),
      panel.background = element_rect(fill = NA, color = "black", linewidth = 1)
    ) +
    guides(fill = guide_legend()) +
    facet_grid(~ Gender, scales = "free", space = "free", labeller = labeller(Gender = c("F" = "Female", "M" = "Male")))


  p
}

peptide_data |>
  imap(~ comparative_reverse_krsa(.x, .y)) |>
  imap(
    ~ ggsave(
      filename = str_glue("{.y}_subset_recombinant_peptide_activity_barplot.png"),
      plot = .x,
      path = "figures/recombinant_peptide",
      width = 12,
      height = 8,
      units = "in",
      bg = "white"
    )
  )
