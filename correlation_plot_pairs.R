# Generate pairwise correlation plot with subject mapping
# Using shared utility functions where possible

# Source shared functions
source("utils/data_processing_functions.R")
library(patchwork)

pair_map_file <- file.path("kinome_data", "subject_pairs.csv")

pair_map_data <- pair_map_file |>
    readr::read_csv(show_col_types = FALSE) |>
    tidyr::drop_na() |>
    dplyr::mutate(
        CTL = dplyr::if_else(CTL == "CTL_832", "CTL_834", CTL),
        value = paste(SCZ, CTL, sep = "_"),
        Pair = stringr::str_remove(Pair, "air")
    ) |>
    dplyr::select(value, Pair)

labels <- tibble::deframe(pair_map_data)

correlation_file <- file.path("results", "combined_score_correlations.csv")

pairwise_data <- correlation_file |>
    readr::read_csv(show_col_types = FALSE) |>
    dplyr::filter(
        stringr::str_detect(var1, "SCZ_\\d{1,3}"),
        stringr::str_detect(var2, "SCZ_\\d{1,3}")
    ) |>
    dplyr::select(X = var1, Y = var2, Correlation = cor) |>
    tidyr::pivot_longer(c(X, Y)) |>
    dplyr::inner_join(pair_map_data, by = "value") |>
    dplyr::rename(Tag = value, value = Pair) |>
    dplyr::select(-Tag) |>
    tidyr::pivot_wider(names_from = name, values_from = value) |>
    tidyr::unnest(cols = c(X, Y))

g <- ggplot2::ggplot(pairwise_data, ggplot2::aes(x = X, y = Y, fill = Correlation))

p <- g +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::scale_x_discrete(name = NULL, labels = labels, expand = c(0L, 0L)) +
    ggplot2::scale_y_discrete(name = NULL, labels = labels, expand = c(0L, 0L)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        text = ggplot2::element_text(size = 18L)
    ) +
    theme(legend.position = "bottom")

ggplot2::ggsave("pairwise_correlations.png",
    width = 10L,
    height = 5L,
    units = "in",
    path = "figures",
    plot = p,
    bg = "white"
)


ggplot2::ggsave("pairwise_correlations.svg",
    width = 18L,
    height = 5L,
    units = "in",
    path = "figures",
    plot = p,
    bg = "transparent"
)
