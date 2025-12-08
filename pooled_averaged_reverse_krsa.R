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
    bind_rows() |>
    left_join(pair_map_data, by = "Group") |>
    select(-Group) |>
    rename(Group = Type)

processed_mapped <- processed |>
    rename(PeptideID = Peptide, Diagnosis = Group) |>
    inner_join(peptide_map, by = "PeptideID", relationship = "many-to-many") |>
    select(-Pair)

get_summary_stats_diff <- function(dataset) {
    dataset |>
        mutate(Group = str_extract(Group, "(SCZ|CTL)")) |>
        pivot_wider(names_from = Group, values_from = Slope) |>
        mutate(Diff = SCZ - CTL) |>
        get_summary_stats(Diff)
}

calculate_pair_differences <- function(peptide_data, kinase) {
    selected_peptides <- peptide_map |>
        filter(Group == !!kinase) |>
        pull(PeptideID) |>
        unique()

    calculated_data <- peptide_data |>
        filter(Peptide %in% selected_peptides) |>
        mutate(Set = "X") |>
        nest(.by = Set) |>
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
            Significant = p < 0.05,
            Normal = normal_p >= 0.01
        ) |>
        select(
            Estimate = estimate, p, conf.low, conf.high, Min = summary_min,
            Max = summary_max, Median = summary_median, Q1 = summary_q1,
            Q3 = summary_q3, SE = summary_se, Effectsize = effect_effsize,
            Magnitude = effect_magnitude, Significant, Normal
        )
}

differed <- c("JNK", "ERK", "P38") |>
    set_names() |>
    map(~ calculate_pair_differences(processed, .x))


g <- processed_mapped |>
    ggplot(aes(x = Diagnosis, y = Slope, fill = Diagnosis))

g + geom_boxplot(width = 0.4) +
    geom_jitter(width = 0.05) +
    scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
    scale_x_discrete(labels = c(CTL = "Control", SCZ = "Schizophrenia")) +
    scale_fill_viridis_d() +
    ylab("Activity") + xlab("") +
    theme_minimal() +
    guides(color = "none") +
    facet_wrap(facets = vars(Group))


ggsave("pooled_reverse_krsa.png", bg = "white", width = 7, height = 7, units = "in")
