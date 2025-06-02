# Develop a library to calculate a combined score for each kinase

library(creedenzymatic)
library(tidyverse)


kinome <- kinome_mp_file

calculate_weights <- function() {
    kinome <- kinome_mp_file

    total_kinases <- nrow(kinome)
    krsa <- kinome |>
        pull(krsa_id) |>
        keep(~ !is.na(.x)) |>
        length()
    uka <- kinome |>
        pull(uka_id) |>
        keep(~ !is.na(.x)) |>
        length()
    kea3 <- kinome |>
        pull(kea3_id) |>
        keep(~ !is.na(.x)) |>
        length()
    ptmsea <- kinome |>
        pull(ptmsea_id) |>
        keep(~ !is.na(.x)) |>
        length()

    coverage <- c(
        KRSA = round(krsa / total_kinases, 6L),
        UKA = round(uka / total_kinases, 6L),
        KEA3 = round(kea3 / total_kinases, 6L),
        PTMSEA = round(ptmsea / total_kinases, 6L)
    )

    normalized_coverage <- round(coverage / sum(coverage), 4L)

    list(
        coverage = coverage,
        normalized_coverage = normalized_coverage
    )
}

penalty_scale <- seq(0L, 1L, 0.05)

method_counts <- c(374L, 377L, 530L, 192L)

penalty_scale_df <- tibble(Scale = penalty_scale) |>
    expand_grid(coverage = calculate_weights()$coverage) |>
    left_join(
        as.data.frame(calculate_weights()) |> rownames_to_column("Method"),
        by = "coverage"
    ) |>
    rename(Coverage = coverage, Weight = normalized_coverage) |>
    select(Method, Coverage, Weight, Scale) |>
    mutate(
        PenaltyC = 1L - (Coverage * Scale),
        PenaltyW = 1L - (Weight * Scale)
    ) |>
    pivot_longer(starts_with("Penalty"), names_to = "PenaltyType", values_to = "Penalty")

g <- ggplot(penalty_scale_df, aes(
    x = Scale, y = Penalty,
    color = Method, group = Method,
    shape = PenaltyType
))

p <- g +
    geom_point() +
    geom_line() +
    theme_bw() +
    facet_wrap(~PenaltyType)
p

penalty_impact_df <- penalty_scale_df |>
    filter(PenaltyType == "PenaltyW") |>
    select(Method, Penalty, Scale) |>
    pivot_wider(names_from = Method, values_from = "Penalty") |>
    mutate(
        b0000 = 1L,
        b0001 = PTMSEA,
        b0010 = KEA3,
        b0011 = KEA3 * PTMSEA,
        b0100 = UKA,
        b0101 = UKA * PTMSEA,
        b0110 = UKA * KEA3,
        b0111 = UKA * KEA3 * PTMSEA,
        b1000 = KRSA,
        b1001 = KRSA * PTMSEA,
        b1010 = KRSA * KEA3,
        b1011 = KRSA * KEA3 * PTMSEA,
        b1100 = KRSA * UKA,
        b1101 = KRSA * UKA * PTMSEA,
        b1110 = KRSA * UKA * KEA3,
        b1111 = 0L
    ) |>
    select(-c(KRSA, UKA, KEA3, PTMSEA)) |>
    pivot_longer(starts_with("b"), names_to = "Missing", values_to = "TotalPenalty")

g2 <- ggplot(penalty_impact_df, aes(x = Missing, y = Scale, fill = TotalPenalty))

p2 <- g2 +
    geom_tile() +
    scale_y_continuous(breaks = seq(0L, 1L, 0.25), minor_breaks = seq(0L, 1L, 0.05)) +
    scale_fill_gradient(high = "red", low = "blue") +
    theme_minimal()

p2

weighted_average_with_coverage <- function(
    p1, p2, p3, p4,
    neutral_value = 0.5,
    penalty_scale = 0.2) {
    calculated_weights_coverage <- calculate_weights()


    calculated_weights <- calculated_weights_coverage$normalized_coverage

    # Compute penalty factors: 1 - (coverage × penalty_scale)
    penalties <- 1L - (calculated_weights_coverage$coverage * penalty_scale)

    # Initialize alpha
    alpha <- 1L

    # Replace missing values (-1) with neutral and apply penalties
    values <- c(p1, p2, p3, p4)
    for (i in 1L:4L) {
        if (values[i] == -1L) {
            alpha <- alpha * penalties[i]
            values[i] <- neutral_value
        }
    }

    # Weighted average
    weighted_avg <- sum(calculated_weights * values)

    # Apply penalty
    final_score <- alpha * weighted_avg
    return(final_score)
}

rescale_combined <- function(x) {
    if (all(x == x[1L])) {
        return(rep(0.5, length(x)))
    } # Handle constant vectors
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

creeden_data <- read_csv(file.path("results", "CTL_F_CTL_M_STK_creedenzymatic.csv")) |>
    mutate(Method = if_else(Method == "PTM-SEA", "PTMSEA", Method)) |>
    select(Kinase, HGNC = hgnc_symbol, Method, Perc) |>
    mutate(Score = Perc) |>
    pivot_wider(names_from = Method, values_from = Score, values_fill = -1L, values_fn = unique) |>
    mutate(
        CombinedScore = pmap_dbl(
            list(KRSA, UKA, KEA3, PTMSEA),
            ~ weighted_average_with_coverage(..1, ..2, ..3, ..4)
        ),
        Percentile = ntile(CombinedScore, 100L),
        Rescaled = rescale_combined(CombinedScore)
    ) |>
    write_csv("test_combined_score.csv")
