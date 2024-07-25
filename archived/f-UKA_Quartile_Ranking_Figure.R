# Green Monster for the KRSA Rankings

library(tidyverse)

colors <- rev(c('#ffffcc', '#c2e699', '#78c679', '#238443'))

core_dataset <- read_csv("results/UKA_quartile_rankings_all.csv") |>
  mutate(quartile = factor(str_c("Q", quartile)))
core_ranks <-
  read_csv("results/UKA_quartile_rankings_aggregated_ordered_all.csv")

excluded <-
  core_dataset |>
  group_by(Kinase) |>
  count() |>
  filter(n < 19) |>
  pull(Kinase)

clean_dataset <- core_dataset |>
  filter(!Kinase %in% excluded)

clean_ranks <- core_ranks |>
  filter(!Kinase %in% excluded)

g <-
  ggplot(clean_dataset, aes(x = dataset, y = Kinase, fill = quartile))

p <- g + geom_tile(size = 0.15, color = "white") +
  scale_fill_manual(
    limits = c("Q1", "Q2", "Q3", "Q4"),
    labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4"),
    values = colors,
    name = ""
  ) +
  theme_minimal() +
  scale_y_discrete(limits = rev(clean_ranks$Kinase)) +
  ggtitle("Quartile Ranking Across the Paired Samples Data for Individual Kinases") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  "UKA_Quartile_Rank_by_Dataset.png",
  plot = p,
  path = "figures",
  width = 15,
  height = 10,
  units = "in",
  bg = "white"
)
