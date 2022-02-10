# Green Monster for the KRSA Rankings

library(tidyverse)

colors <- rev(c('#ffffcc', '#c2e699', '#78c679', '#238443'))

core_dataset <- read_csv("results/KRSA_quartile_rankings_all.csv") |>
  mutate(quartile = factor(str_c("Q", quartile))) |>
  filter(Kinase != "PDK")
core_ranks <- read_csv("results/KRSA_quartile_rankings_aggregated_ordered_all.csv") |>
  filter(Kinase != "PDK")

g <-
  ggplot(core_dataset, aes(x = dataset, y = Kinase, fill = quartile))

p <- g + geom_tile(size = 0.15, color = "white") +
  scale_fill_manual(
    limits = c("Q1", "Q2", "Q3", "Q4"),
    labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4"),
    values = colors,
    name = ""
  ) +
  theme_minimal() +
  scale_y_discrete(limits = rev(core_ranks$Kinase)) +
  ggtitle("Quartile Ranking Across the Paired Samples Data for Kinase Families") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("KRSA_Quartile_Rank_by_Dataset.png", plot = p, path = "figures", width = 15, height = 10, units = "in", bg = "white")
