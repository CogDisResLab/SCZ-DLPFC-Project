# Generate reverse KRSA Plots for our top 3 families

suppressPackageStartupMessages({
    library(tidyverse)
    library(rstatix)
    library(patchwork)
})

signal_files <- list.files("results", "cell.*signal") |>
    discard(~ str_detect(.x, fixed("pooled")))

peptide_map <- read_csv(file.path("reference_data", "MAPK-mapped-kinases.csv"))

pair_map_file <- file.path("kinome_data", "subject_pairs.csv")

pair_map_data <- pair_map_file |>
    readr::read_csv(show_col_types = FALSE) |>
    tidyr::drop_na() |>
    dplyr::mutate(
        CTL = dplyr::if_else(CTL == "CTL_832", "CTL_834", CTL),
        Pair = stringr::str_remove(Pair, fixed("air"))
    ) |>
    dplyr::select(CTL, SCZ, Pair) |>
    pivot_longer(c(CTL, SCZ), names_to = "Type", values_to = "Group")

process_file <- function(filename) {
    filepath <- file.path("results", filename)

    signal_data <- read_csv(filepath) |>
        filter(Peptide %in% peptide_map$PeptideID) |>
        select(Group, Peptide, Slope = slope)
}

processed <- signal_files |>
    process_file() |>
    bind_rows()

get_summary_stats_diff <- function(dataset) {
    dataset |>
        mutate(Group = str_extract(Group, "(SCZ|CTL)")) |>
        pivot_wider(names_from = Group, values_from = Slope) |>
        mutate(Diff = SCZ - CTL) |>
        get_summary_stats(Diff)
}

calculate_pair_differences <- function(peptide_data, group) {
    selected_peptides <- peptide_map |>
        filter(Group == !!group) |>
        pull(PeptideID) |>
        unique()

    calculated_data <- peptide_data |>
        left_join(pair_map_data, by = "Group") |>
        select(-Type) |>
        filter(Peptide %in% selected_peptides) |>
        mutate(Group = as_factor(str_extract(Group, "(SCZ|CTL)"))) |>
        nest(.by = "Pair") |>
        mutate(
            normal = map(data, ~ shapiro_test(.x, Slope)),
            # NOTE: rstatix::wilcox_test reports (ref.group - comparison).
            # We want SCZ - CTL, so ref.group must be "SCZ".
            test = map(
                data, ~ wilcox_test(Slope ~ Group,
                    data = .x, paired = TRUE, detailed = TRUE,
                    ref.group = "SCZ"
                )
            ),
            summary = map(data, ~ get_summary_stats_diff(.x)),
            effect = map(data, ~ wilcox_effsize(
                Slope ~ Group,
                data = .x, paired = TRUE,
                ci = TRUE, conf.level = 0.1,
                nboot = 1000L
            ))
        ) |>
        unnest_wider(test) |>
        unnest_wider(normal, names_sep = "_") |>
        unnest_wider(summary, names_sep = "_") |>
        unnest_wider(effect, names_sep = "_") |>
        mutate(
            Significant = p < 0.1,
            Normal = normal_p >= 0.1
        ) |>
        select(
            Pair,
            Estimate = estimate, p, conf.low, conf.high, Min = summary_min,
            Max = summary_max, Median = summary_median, Q1 = summary_q1,
            Q3 = summary_q3, SE = summary_se, Effectsize = effect_effsize,
            Magnitude = effect_magnitude, Significant, Normal
        )
}

differed <- c("JNK", "ERK", "P38") |>
    set_names() |>
    map(~ calculate_pair_differences(processed, .x))

comparative_reverse_krsa <- function(kinase, dataset) {
    plot_data <-
        mutate(dataset, Lower = Estimate - SE, Upper = Estimate + SE)

    g <- ggplot(
        plot_data,
        aes(
            x = Pair, y = Estimate,
            lower = Q1, upper = Q3,
            middle = Estimate, ymin = conf.low, ymax = conf.high,
            fill = Significant, color = Significant
        )
    )

    order <- plot_data |>
        arrange(desc(Effectsize)) |>
        pull(Pair)

    p <- g +
        geom_hline(yintercept = 0L, color = "grey50", lwd = 0.5) +
        geom_hline(yintercept = 0.2, color = "grey80", lwd = 1L) +
        geom_hline(yintercept = -0.2, color = "grey80", lwd = 1L) +
        geom_point(aes(size = Effectsize)) +
        geom_errorbar(width = 0.3) +
        theme_minimal() +
        ggtitle(str_glue("Comparison of datasets for {kinase} family")) +
        scale_color_manual(limits = c(TRUE, FALSE), values = c("red", "black")) +
        scale_fill_manual(limits = c(TRUE, FALSE), values = c("red", "black")) +
        scale_y_continuous(
            name = expression(Log[2L] ~ Fold ~ Change),
            limits = c(-1.5, 1.5),
            breaks = seq(-2.6, 2.6, 0.2)
        ) +
        scale_x_discrete(name = "Subject Pair", limits = order) +
        scale_size_continuous(limits = c(0L, 1.1), name = "Effect Size", breaks = seq(0L, 1L, 0.2)) +
        theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24L)) +
        guides(fill = "none", color = "none")


    p
}

plots <- differed |>
    imap(~ comparative_reverse_krsa(.y, .x))

plots |>
    imap(~ ggsave(str_glue("{.y}_LFC_Errorbar_Plot.svg"),
        plot = .x, width = 10L,
        height = 7L, units = "in", path = "figures"
    ))

combined_plot <- (plots$P38 / plots$JNK / plots$ERK) +  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("combined_errorbar_plot.png", plot = combined_plot, height = 24, width = 18L, units = "in", path = "figures")
ggsave("combined_errorbar_plot.svg", plot = combined_plot, height = 24, width = 18L, units = "in", path = "figures", bg = "transparent")
