# Green Monster for the KRSA Rankings

library(tidyverse)

colors <- rev(c("#ffffcc", "#c2e699", "#78c679", "#238443"))

calculate_label_pair <- function(name) {
  run_number <- str_extract(name, "R\\d") |>
    str_remove("R") |>
    as.integer()
  pair_number <- str_extract(name, "P\\d") |>
    str_remove("P") |>
    as.integer()

  number <- (run_number - 1L) * 6L + pair_number

  number |> as.character() |> str_pad(2L, pad = "0")
}


core_dataset <- list.files("results", pattern = "krsa_table_R\\dC\\dP\\d.+.csv", full.names = TRUE) |>
  set_names(~ basename(.x) |> str_extract("(R\\dC\\dP\\d_.+)_STK.csv", 1L)) |>
  map(~ read.csv(.x)) |>
  map(~ select(.x, Kinase, Score = AvgZ)) |>
  map(~ distinct(.x)) |>
  map(~ arrange(.x, desc(abs(Score)))) |>
  map(~ mutate(
    .x,
    Rank = row_number(desc(abs(Score))),
    Quartile = ntile(desc(abs(Score)), 4L)
  )) |>
  bind_rows(.id = "dataset") |>
  mutate(
    Label = str_extract(dataset, "R\\dC\\dP\\d"),
    LabelPair = calculate_label_pair(Label),
    Quartile_factor = factor(str_c("Q", Quartile))
  ) |>
  select(Label, LabelPair, Kinase, Score, Quartile, Quartile_factor)

ranked_data <- core_dataset |>
  select(Label, Kinase, Quartile) |>
  distinct() |>
  pivot_wider(names_from = Label, values_from = Quartile, values_fill = 100L) |>
  pivot_longer(cols = -Kinase, names_to = "Label", values_to = "Quartile") |>
  summarise(MeanQuartile = mean(Quartile, na.rm = TRUE), .by = Kinase) |>
  arrange(MeanQuartile)

g <- ggplot(core_dataset, aes(x = LabelPair, y = Kinase, fill = Quartile_factor))

p <- g + geom_tile() +
  scale_fill_manual(
    limits = c("Q1", "Q2", "Q3", "Q4"),
    labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4"),
    values = colors,
    name = ""
  ) +
  theme_minimal() +
  scale_y_discrete(limits = rev(ranked_data$Kinase)) +
  guides(fill = "none")

ggsave("KRSA_Quartile_Rank_by_Dataset.png", plot = p, path = "figures", width = 6.5, height = 9L, units = "in", bg = "white")

ggsave("KRSA_Quartile_Rank_by_Dataset.svg", plot = p, path = "figures", width = 6.5, height = 9L, units = "in", bg = "white")
