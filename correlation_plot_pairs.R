library(tidyverse)

pair_map_file <- file.path("kinome_data", "subject_pairs.csv")

pair_map_data <- pair_map_file |>
    read_csv() |>
    na.omit() |>
    mutate(
        CTL = if_else(CTL == "CTL_832", "CTL_834", CTL),
        value = paste(SCZ, CTL, sep = "_"),
        Pair = str_remove(Pair, "air")
    ) |>
    select(value, Pair)

labels <- deframe(pair_map_data)

correlation_file <- file.path("results", "combined_score_correlations.csv")


pairwise_data <- correlation_file |>
    read_csv() |>
    filter(
        str_detect(var1, "SCZ_\\d{1,3}"), # nolint: nonportable_path_linter.
        str_detect(var2, "SCZ_\\d{1,3}") # nolint: nonportable_path_linter.
    ) |>
    select(X = var1, Y = var2, Correlation = cor) |>
    pivot_longer(c(X, Y)) |>
    inner_join(pair_map_data) |>
    rename(Tag = value, value = Pair) |>
    select(-Tag) |>
    pivot_wider(names_from = name, values_from = value) |>
    unnest(cols = c(X, Y))

g <- ggplot(pairwise_data, aes(x = X, y = Y, fill = Correlation))

p <- g + geom_tile() +
    scale_fill_viridis_c() +
    scale_x_discrete(name = NULL, labels = labels, expand = c(0L, 0L)) +
    scale_y_discrete(name = NULL, labels = labels, expand = c(0L, 0L)) +
    theme_minimal() +
    theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), text = element_text(size = 18L))

ggsave("pairwise_correlations.png",
    width = 10L,
    height = 5L,
    units = "in",
    path = "figures", plot = p, bg = "white"
)
