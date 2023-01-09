# Visualize the peptide-level data

library(tidyverse)

files <-
  list.files(pattern = "*filtered_peptides_data.csv", recursive = TRUE)

kinases <- str_extract(files, "\\w+_filtered") |>
  str_remove("_filtered") |>
  set_names()

peptide_data <- files |>
  set_names(kinases) |>
  map(~ read_csv(.x))

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

peptide_data |>
  map(generate_heatmap) |>
  map2(
    kinases,
    ~ ggsave(
      filename = str_glue("{.y}_recombinant_peptide_activity.png"),
      plot = .x,
      path = "figures/recombinant_peptide",
      width = 12,
      height = 8,
      units = "in",
      bg = "white"
    )
  )
