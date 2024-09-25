# Green Monster for the KRSA Rankings

library(tidyverse)

colors <- rev(c("#ffffcc", "#c2e699", "#78c679", "#238443"))

core_dataset <- list.files("results", pattern = "krsa_table_R\\dC\\dP\\d.+.csv", full.names = TRUE) |>
  set_names(~ basename(.x) |> str_extract("(R\\dC\\dP\\d_.+)_STK.csv", 1)) |>
  map(~ read.csv(.x)) |>
  map(~ select(.x, Kinase, Score = AvgZ)) |>
  map(~ distinct(.x)) |>
  map(~ arrange(.x, desc(abs(Score)))) |>
  map(~ mutate(
    .x,
    Rank = row_number(desc(abs(Score))),
    Quartile = ntile(desc(abs(Score)), 4)
  )) |>
  bind_rows(.id = "dataset") |>
  mutate(
    Label = str_extract(dataset, "R\\dC\\dP\\d"),
    Quartile_factor = factor(str_c("Q", Quartile))
  ) |>
  select(Label, Kinase, Score, Quartile, Quartile_factor)

ranked_data <- core_dataset |>
  select(Label, Kinase, Quartile) |>
  distinct() |>
  pivot_wider(names_from = Label, values_from = Quartile, values_fill = 100) |>
  pivot_longer(cols = -Kinase, names_to = "Label", values_to = "Quartile") |>
  summarise(MeanQuartile = mean(Quartile, na.rm = TRUE), .by = Kinase) |>
  arrange(MeanQuartile)

g <-  ggplot(core_dataset, aes(x = Label, y = Kinase, fill = Quartile_factor))

p <- g + geom_raster() +
  scale_fill_manual(
    limits = c("Q1", "Q2", "Q3", "Q4"),
    labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4"),
    values = colors,
    name = ""
  ) +
  theme_minimal() +
  scale_y_discrete(limits = rev(ranked_data$Kinase)) +
  ggtitle("Quartile Ranking Across the Paired Samples Data for Kinase Families") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("KRSA_Quartile_Rank_by_Dataset.png", plot = p, path = "figures", width = 6.5, height = 9, units = "in", bg = "white")

ggsave("KRSA_Quartile_Rank_by_Dataset.svg", plot = p, path = "figures", width = 6.5, height = 9, units = "in", bg = "white")
