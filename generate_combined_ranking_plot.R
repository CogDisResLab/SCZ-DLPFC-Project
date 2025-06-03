# Green Monster for the Combined Rankings

library(tidyverse)

colors <- rev(c("#ffffcc", "#c2e699", "#78c679", "#238443"))
pair_df <- read_csv(file.path("kinome_data", "subject_pairs.csv"))

calculate_label_pair <- function(name) {
  scz <- str_extract(name, "SCZ_\\d+")
  ctl <- str_extract(name, "CTL_\\d+")

  pair_df |>
    filter(SCZ == scz, CTL == ctl) |>
    summarise(Pair = first(Pair, default = NA_character_)) |>
    pull(Pair)
}


core_dataset <- list.files("results",
  pattern = "creedencombined", full.names = TRUE
) |>
  keep(~ str_detect(.x, "SCZ_\\d+_CTL_\\d+")) |> # nolint: nonportable_path_linter.
  set_names(~ basename(.x) |> str_remove("_STK_creedencombined.csv")) |>
  map(~ read.csv(.x)) |>
  map(~ select(.x, HGNC, Score = CombinedScore)) |>
  map(~ distinct(.x)) |>
  map(~ arrange(.x, desc(abs(Score)))) |>
  map(~ mutate(
    .x,
    Rank = row_number(desc(abs(Score))),
    Quartile = ntile(desc(abs(Score)), 4L)
  )) |>
  bind_rows(.id = "dataset") |>
  mutate(
    Label = dataset,
    LabelPair = Label,
    Quartile_factor = factor(str_c("Q", Quartile))
  ) |>
  select(Label, LabelPair, HGNC, Score, Quartile, Quartile_factor)

ranked_data <- core_dataset |>
  select(Label, HGNC, Quartile) |>
  distinct() |>
  pivot_wider(names_from = Label, values_from = Quartile, values_fill = 100L) |>
  pivot_longer(cols = -HGNC, names_to = "Label", values_to = "Quartile") |>
  summarise(MeanQuartile = mean(Quartile, na.rm = TRUE), .by = HGNC) |>
  arrange(MeanQuartile)

g <- ggplot(core_dataset, aes(x = LabelPair, y = HGNC, fill = Score))

p <- g + geom_tile() +
  scale_fill_gradient(low = "yellow", high = "red") +
  theme_minimal() +
  scale_y_discrete(limits = rev(ranked_data$HGNC)) +
  guides(fill = "none")

ggsave("Combined_Quartile_Rank_by_Dataset.png",
  plot = p, path = "figures",
  width = 6.5, height = 9L * 5, units = "in", bg = "white"
)

ggsave("Combined_Quartile_Rank_by_Dataset.svg",
  plot = p, path = "figures",
  width = 6.5, height = 9L * 5, units = "in", bg = "white"
)
