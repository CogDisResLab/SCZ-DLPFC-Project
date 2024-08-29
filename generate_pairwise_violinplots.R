# Create pairwise heatmaps for each matched pair of samples
#
# We will scale globally so all heatmaps have the same scale

library(tidyverse)
library(ggpubr)
library(patchwork)

signal_files <- list.files("results", "cell.*signal", full.names = TRUE) |>
  discard(~ str_detect(.x, fixed("SCZ_11"))) |>
  set_names(~ str_extract(basename(.x), "SCZ_\\d+_CTL_\\d+")) # nolint: nonportable_path_linter.

pair_numbers <- read_csv("kinome_data/subject_pairs.csv") |>
  select(Pair, SCZ, CTL) |>
  mutate(
    CTL = str_replace(CTL, fixed("832"), fixed("834")),
    ID = str_c(SCZ, CTL, sep = "_")
  ) |>
  select(-SCZ, -CTL, ID, Pair) |>
  filter(!is.na(Pair)) |>
  select(ID, Pair) |>
  deframe()

signal_files <- set_names(signal_files, pair_numbers)

col_spec <- cols(
  Group = col_character(),
  Barcode = col_skip(),
  SampleName = col_skip(),
  Peptide = col_character(),
  slope = col_double(),
  r.seq = col_skip(),
  nobs = col_skip()
)

create_violinplot <- function(data, comparison) {
  clean_data <- data |>
    mutate(
      Group = str_extract(Group, "SCZ|CTL")
    )

  ggviolin(
    clean_data,
    "Group",
    "slope",
    xlab = NULL,
    ylab = NULL,
    fill = "Group",
    palette = "uchicago",
    ylim = scale_limits,
    legend = "none",
    yticks.by = 2L,
    add = "boxplot",
    add.params = list(
      fill = "white",
      color = "black"
    )
  ) +
    geom_bracket(
      aes(
        label = p.signif,
        xmin = group1,
        xmax = group2
      ),
      data = compare_means(
        slope ~ Group,
        data = clean_data, method = "t.test"
      ),
      y.position = max(scale_limits) * 0.95
    ) +
    labs(caption = comparison, x = NULL, y = NULL) +
    theme(
      plot.caption = element_text(hjust = 0.5)
    )
}

signal_data <- signal_files |>
  map(~ read_csv(.x, col_types = col_spec)) |>
  bind_rows(.id = "Comparison")

scale_limits <- c(
  floor(min(signal_data[["slope"]])),
  ceiling(max(signal_data[["slope"]]))
)

breaks <- seq(
  scale_limits[[1L]],
  scale_limits[[2L]],
  length.out = 6L
)

violinplots <- signal_data |>
  nest(.by = Comparison) |>
  mutate(
    violinplot = map2(
      data,
      Comparison,
      ~ create_violinplot(.x, .y)
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
    width = 5L,
    height = 6.68,
    path = "figures",
    dpi = 600L,
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
    width = 5L,
    height = 6.68,
    path = "figures",
    dpi = 600L,
    bg = "white",
    create.dir = TRUE
  ))

combined_plot <- wrap_plots(violinplots[-c(10L)], nrow = 3L, ncol = 6L) # nolint: unnecessary_concatenation_linter.

ggsave(
  filename = file.path("figures", "combined_violinplots.png"),
  plot = combined_plot,
  width = 7.5,
  height = 5L,
  dpi = 300L,
  bg = "white"
)

ggsave(
  filename = file.path("figures", "combined_violinplots.svg"),
  plot = combined_plot,
  width = 7.5,
  height = 5L,
  dpi = 300L,
  bg = "white"
)
