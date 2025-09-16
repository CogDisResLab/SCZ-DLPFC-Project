# Green Monster for the KRSA Rankings

library(tidyverse)

colors <- rev(c("#ffffcc", "#c2e699", "#78c679", "#238443"))

signal_files <- list.files("results", "cell.*krsa", full.names = TRUE) |>
    discard(~ str_detect(.x, fixed("SCZ_11"))) |>
    discard(~ str_detect(.x, fixed("pooled"))) |>
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

core_dataset <- signal_files |>
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
        Label = {
            dataset |> str_remove("air")
        },
        LabelPair = dataset,
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

g <- ggplot(core_dataset, aes(x = Label, y = Kinase, fill = Quartile_factor))

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

ggsave("KRSA_Quartile_Rank_by_Dataset.svg", plot = p, path = "figures", width = 6.5, height = 9L, units = "in", bg = "transparent")
