# Reverse KRSA on the rats data

library(tidyverse)
library(KRSA)
library(gt)
library(knitr)

signal_file <- "raw/rats-haloperidol/rats-haloperidol_SigmBg.txt"
saturation_file <- "raw/rats-haloperidol/rats-haloperidol_Saturation.txt"

dataset_max <- krsa_read(signal_file, saturation_file) |>
  filter(SampleName == "Brain") |>
  krsa_qc_steps() |>
  mutate(Group = if_else(`Comment 1` == "CT", "CTL", "HPD")) |>
  krsa_extractEndPointMaxExp(type = "STK")

dataset_all <- krsa_read(signal_file, saturation_file) |>
  filter(SampleName == "Brain") |>
  krsa_qc_steps() |>
  mutate(Group = if_else(`Comment 1` == "CT", "CTL", "HPD")) |>
  krsa_extractEndPoint(type = "STK")

scaled_model <- dataset_max |>
  krsa_filter_lowPeps(threshold = 5) |>
  {
    \(x) krsa_scaleModel(dataset_all, x)
  }()

process_data <- function(signal_file, saturation_file) {
  dataset_max <- krsa_read(signal_file, saturation_file) |>
    filter(str_detect(SampleName, "CTL") |
             str_detect(SampleName, "SCZ")) |>
    krsa_qc_steps() |>
    mutate(Group = str_extract(SampleName, "^\\w{3}")) |>
    krsa_extractEndPointMaxExp(type = "STK")

  dataset_all <- krsa_read(signal_file, saturation_file) |>
    filter(str_detect(SampleName, "CTL") |
             str_detect(SampleName, "SCZ")) |>
    krsa_qc_steps() |>
    mutate(Group = str_extract(SampleName, "^\\w{3}")) |>
    krsa_extractEndPoint(type = "STK")

  scaled_model <- dataset_max |>
    krsa_filter_lowPeps(threshold = 5) |>
    {
      \(x) krsa_scaleModel(dataset_all, x)
    }()

  out <- list(
    dataset_max_exp = dataset_max,
    dataset_all_exp = dataset_all,
    scaled_models = scaled_model$scaled
  )
}

generate_group_differences <-
  function(ctl, scz, dataset_max_exp, scaled_models) {
    comparisons <- list(Comp1 = c("SCZ", "CTL"))
    sample_ids <- c(ctl, scz)

    selected_peptides_compared <-
      krsa_quick_filter(
        data = dataset_max_exp,
        data2 = scaled_models,
        signal_threshold = 5,
        r2_threshold = 0.7,
        samples = sample_ids,
        groups = comparisons$Comp1
      )

    group_diff_df <-
      krsa_group_diff(
        scaled_models,
        comparisons$Comp1,
        samples = sample_ids,
        selected_peptides_compared,
        byChip = T
      )

    group_diff_df
  }


generate_reverse_krsa_ggplot <-
  function(group_diff, hits) {
    fig <-
      krsa_reverse_krsa_plot(KRSA_coverage_STK_PamChip_87102_v2,
                             group_diff,
                             hits,
                             0.2,
                             F)

    fig

  }
