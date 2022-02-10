# Create a CMGC Volcano Plot

library(tidyverse)

dataset <- read_csv("results/CMGC_subset_SCZ_lookup.csv") |>
  mutate(indicator = case_when(
    P_Value < 0.05 & Log2FC > 1 ~ "HiSig",
    P_Value < 0.05 ~ "Sig",
    Log2FC > 1 ~ "Hi",
    TRUE ~ "None"
  ),
  P_Value = if_else(P_Value == 0, 0.000002, P_Value))


g <- ggplot(dataset, aes(x = Log2FC, y = -log10(P_Value), color = indicator))

p <- g + geom_point() +
  theme_minimal() +
  scale_color_manual(
    breaks = c("HiSig", "Sig", "Hi", "None"),
    labels = c("Overexpressed and Significant", "Signficant Only", "Overexpressed Only", "Unremarkable"),
    values = c("red", "blue", "black", "grey"),
    name = "Status"
  ) +
  scale_x_continuous(breaks = seq(-3, 3, 1)) +
  scale_y_continuous(breaks = seq(0, 6, 1),
                     limits = c(0, 6)) +
  theme(legend.position = "bottom")

ggsave("CMGC-Volcano-Plot.png", plot = p, width = 16, height = 8, units = "in", path = "figures", bg = "white")

