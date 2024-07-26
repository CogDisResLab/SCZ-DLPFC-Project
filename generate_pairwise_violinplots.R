# Create pairwise heatmaps for each matched pair of samples
#
# We will scale globally so all heatmaps have the same scale

library(tidyverse)

signal_files <- list.files("results", "cell.*signal", full.names = TRUE) |>
  discard(~ str_detect(.x, fixed("SCZ_11"))) |>
  set_names(~ str_extract(basename(.x), "SCZ_\\d+_CTL_\\d+")) # nolint: nonportable_path_linter.

col_spec <- cols(
  Group = col_character(),
  Barcode = col_skip(),
  SampleName = col_skip(),
  Peptide = col_character(),
  slope = col_double(),
  r.seq = col_skip(),
  nobs = col_skip()
)

create_violinplot <- function(data) {
  ggplot(data, aes(x = Group, y = slope, fill = Group)) +
    geom_violin() +
    geom_jitter(width = 0.1) +
    scale_x_discrete(labels = \(x) str_extract(x, "SCZ|CTL")) +
    scale_fill_discrete(type = hcl.colors(2L, palette = "Harmonic")) +
    xlab("Phenotype") +
    ylab("Signal Intensity") +
    guides(fill = "none") +
    ggprism::theme_prism()
}

signal_data <- signal_files |>
  map(~ read_csv(.x, col_types = col_spec)) |>
  bind_rows(.id = "Comparison")

violinplots <- signal_data |>
  nest(.by = Comparison) |>
  mutate(
    violinplot = map(
      data,
      ~ create_violinplot(.x)
    )
  ) |>
  select(Comparison, violinplot) |>
  deframe()

figures_png <- violinplots |>
  imap(~ ggsave(
    filename = file.path(
      "pairwise_violinplots",
      "png",
      str_glue("{.y}_violinplots.png")
    ),
    plot = .x,
    width = 3L,
    height = 5L,
    path = "figures",
    dpi = 300L,
    bg = "white",
    create.dir = TRUE
  ))

figures_svg <- violinplots |>
  imap(~ ggsave(
    filename = file.path(
      "pairwise_violinplots",
      "svg",
      str_glue("{.y}_violinplots.svg")
    ),
    plot = .x,
    width = 3L,
    height = 5L,
    path = "figures",
    dpi = 300L,
    bg = "white",
    create.dir = TRUE
  ))
