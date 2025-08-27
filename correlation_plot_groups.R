library(tidyverse)

labels <- c(
    "C3N_C3NPC", "C3N_D3N", "C3N_D3NPC", "CTL_F_CTL_M", "CTL_SCZ",
    Haloperidol_Control = "Rats", "SCZ_CTL", SCZ_F_CTL_F = "SCZ_F", SCZ_F_SCZ_M",
    "SCZ_M_CTL_M"
)

correlation_file <- file.path("results", "combined_score_correlations.csv")


correlation_data <- correlation_file |>
    read_csv() |>
    filter(
        str_detect(var1, "SCZ_\\d{1,3}", negate = TRUE), # nolint: nonportable_path_linter.
        str_detect(var2, "SCZ_\\d{1,3}", negate = TRUE) # nolint: nonportable_path_linter.
    ) |>
    select(X = var1, Y = var2, Correlation = cor)

g <- ggplot(correlation_data, aes(x = X, y = Y, fill = Correlation))

p <- g + geom_tile() +
    scale_fill_viridis_c() +
    scale_x_discrete(name = NULL, labels = labels, expand = c(0L, 0L)) +
    scale_y_discrete(name = NULL, labels = labels, expand = c(0L, 0L)) +
    theme_minimal() +
    theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), text = element_text(size = 18L))

ggsave("groupwise_correlations.png",
    width = 10L,
    height = 5L,
    units = "in",
    path = "figures", plot = p, bg = "white"
)
