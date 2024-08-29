# Create pairwise heatmaps for each matched pair of samples
#
# We will scale globally so all heatmaps have the same scale

library(tidyverse)
library(ggheatmap)
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

signal_data <- signal_files |>
  map(~ read_csv(.x, col_types = col_spec)) |>
  bind_rows(.id = "Comparison")

comparison_annotation <- signal_data |>
  select(Comparison, Group) |>
  distinct()

signal_data_wide <- signal_data |>
  select(Group, Peptide, slope) |>
  pivot_wider(names_from = Peptide, values_from = slope, values_fill = 0L) |>
  column_to_rownames("Group")

scaled_signal <- signal_data_wide |>
  as.matrix() |>
  scale() |>
  round(5L) |>
  as.data.frame() |>
  rownames_to_column("Group") |>
  pivot_longer(-Group, names_to = "Peptide", values_to = "ScaledValue") |>
  left_join(comparison_annotation, by = "Group") |>
  write_csv(file.path("results", "scaled_signal.csv"))

scale_limits <- c(
  floor(min(scaled_signal[["ScaledValue"]])),
  ceiling(max(scaled_signal[["ScaledValue"]]))
)

heatmaps <- scaled_signal |>
  nest(.by = Comparison) |>
  mutate(
    widened = map(
      data,
      ~ pivot_wider(.x,
        names_from = Group,
        values_from = ScaledValue
      ) |>
        column_to_rownames("Peptide") |>
        as.matrix()
    ),
    heatmap = map2(
      widened,
      Comparison,
      ~ ggheatmap(.x,
        color = hcl.colors(100L, "Red-Green"),
        legendName = "Activity",
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        show_cluster_rows = FALSE,
        text_show_rows = FALSE,
        scale = "none"
      ) +
        scale_x_discrete(labels = \(x) str_extract(x, "SCZ|CTL")) +
        scale_fill_gradientn(
          colours = hcl.colors(100L, "Red-Green"),
          limits = scale_limits
        ) +
        guides(fill = "none") +
        labs(caption = .y) +
        theme(plot.caption = element_text(hjust = 0.5))
    )
  ) |>
  select(Comparison, heatmap) |>
  deframe()

figures_png <- heatmaps |>
  imap(~ ggsave(
    filename = file.path(
      "pairwise_heatmaps",
      "png",
      str_glue("{.y}_heatmaps.png")
    ),
    plot = .x,
    width = 1.25,
    height = 1.67,
    path = "figures",
    dpi = 300L,
    bg = "white",
    create.dir = TRUE
  ))

figures_svg <- heatmaps |>
  imap(~ ggsave(
    filename = file.path(
      "pairwise_heatmaps",
      "svg",
      str_glue("{.y}_heatmaps.svg")
    ),
    plot = .x,
    width = 1.25,
    height = 1.67,
    path = "figures",
    dpi = 300L,
    bg = NULL,
    create.dir = TRUE
  ))

combined_plot <- wrap_plots(heatmaps[-c(10L)], nrow = 3L, ncol = 6L) # nolint: unnecessary_concatenation_linter.

ggsave(
  filename = file.path("figures", "combined_heatmaps.png"),
  plot = combined_plot,
  width = 7.5,
  height = 5L,
  dpi = 300L,
  bg = "white"
)

ggsave(
  filename = file.path("figures", "combined_heatmaps.svg"),
  plot = combined_plot,
  width = 7.5,
  height = 5L,
  dpi = 300L,
  bg = "white"
)
