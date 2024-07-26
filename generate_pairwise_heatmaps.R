# Create pairwise heatmaps for each matched pair of samples
#
# We will scale globally so all heatmaps have the same scale

library(tidyverse)
library(ggheatmap)

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
    heatmap = map(
      widened,
      ~ ggheatmap(.x,
        color = hcl.colors(100L, "Red-Green"),
        legendName = "Activity",
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        show_cluster_rows = FALSE,
        text_show_rows = FALSE,
        scale = "none"
      )
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
    width = 3L,
    height = 5L,
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
    width = 3L,
    height = 5L,
    path = "figures",
    dpi = 300L,
    bg = "white",
    create.dir = TRUE
  ))
