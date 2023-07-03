# Visualize the peptide-level data for specific pairs as barplots

library(tidyverse)
library(ggprism)

files <-
  list.files(pattern = "*filtered_peptides_data.csv", recursive = TRUE)

kinases <- str_extract(files, "\\w+_filtered") |>
  str_remove("_filtered") |>
  set_names()

selected_samples <-
  c("R1C2P3",
    "R2C2P4",
    "R3C3P5",
    "R4C1P2",
    "R1C1P2",
    "R2C2P3",
    "R3C1P1",
    "R3C2P3")

sample_metadata <- read_csv("annotation/sample_matching.csv") |>
  select(Pair = Designation, Gender) |>
  filter(Pair %in% selected_samples)

peptide_data <- files |>
  set_names(kinases) |>
  map( ~ read_csv(.x, show_col_types = FALSE)) |>
  map( ~ inner_join(.x, sample_metadata, by = "Pair"))

comparative_reverse_krsa <- function(dataset, kinase) {
  summarised_data <- dataset |>
    group_by(Pair, Gender) |>
    summarise(
      Mean = mean(Score, na.rm = TRUE),
      SD = sd(Score, na.rm = TRUE),
      N = n(),
      SEM = SD / N,
      .groups = "drop"
    ) |>
    mutate(YMAX = if_else(Mean > 0, Mean + SEM, Mean - SEM))

  sample_order <- summarised_data |>
    ungroup() |>
    arrange(Gender, Pair) |>
    mutate(
      PairNum = row_number(),
      PairLabel = str_glue("Pair{str_pad(PairNum, width=2, pad = 0)}")
    ) |>
    select(PairLabel, Pair) |>
    deframe()

  g <- ggplot(summarised_data,
              aes(x = Pair, y = Mean, fill = Gender))

  p <- g +
    geom_errorbar(aes(x = Pair,
                      ymin = YMAX,
                      ymax = YMAX),
                  width = 0.5,
                  position = "dodge",
                  show.legend = FALSE) +
    geom_linerange(aes(x = Pair,
                       ymin = Mean,
                       ymax = YMAX),
                   show.legend = FALSE) +
    geom_bar(key_glyph = "rect", stat = "identity", show.legend = FALSE) +
    theme_prism() +
    scale_fill_prism(palette = "black_and_white") +
    scale_x_discrete(limits = sample_order, labels = names(sample_order)) +
    scale_y_continuous(
      name = expression(Log[2]~Fold~change),
      limits = c(-0.25, 1.25),
      breaks = seq(-2.5, 2.5, 0.25)
    ) +
    geom_segment(x = 0.75, y = 1.25, xend = 4.25, yend = 1.25, linewidth = 1.25) +
    geom_segment(x = 4.75, y = 1.25, xend = 8.25, yend = 1.25, linewidth = 1.25) +
    geom_text(x = 2.5, y = 1.3, label = "Female", size = 8) +
    geom_text(x = 6.5, y = 1.3, label = "Male", size = 8) +
    xlab("") +
    theme(axis.title.y = element_text(face = "bold"))


  p
}

peptide_data |>
  imap( ~ comparative_reverse_krsa(.x, .y)) |>
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

