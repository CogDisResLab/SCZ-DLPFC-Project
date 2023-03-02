# Generate reverse KRSA Plots for our top 3 families

library(tidyverse)
library(KRSA)
library(gt)
library(knitr)

signal_files <-
  c(
    "raw/Paired-Raw-Data/run1/_KRSAsfiles_Median_SigmBg_210518210643.txt",
    "raw/Paired-Raw-Data/run2/_KRSA_files_Median_SigmBg_210520201325.txt",
    "raw/Paired-Raw-Data/run3/_KRSA_files_Median_SigmBg_210519103637.txt",
    "raw/Paired-Raw-Data/run4/_KRSA_files_Median_SigmBg_210605104306.txt"
  )
saturation_files <-
  c(
    "raw/Paired-Raw-Data/run1/_KRSAsfiles_Signal_Saturation_210518210643.txt",
    "raw/Paired-Raw-Data/run2/_KRSA_files_Signal_Saturation_210520201326.txt",
    "raw/Paired-Raw-Data/run3/_KRSA_files_Signal_Saturation_210519103637.txt",
    "raw/Paired-Raw-Data/run4/_KRSA_files_Signal_Saturation_210605104307.txt"
  )


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



processed <- map2(signal_files, saturation_files,
                  ~ process_data(.x, .y)) |>
  pmap(bind_rows) |>
  map( ~ filter(.x, SampleName != "CTL_859")) |>
  map(~ mutate(.x, SampleName = if_else(
    SampleName == "CTL_834", "CTL_832", SampleName
  )))



samples_ids <- read_csv("annotation/subject_pairs.csv") |>
  select(-Run) |>
  filter(CTL != "CTL_859") |>
  mutate(
    comp_name = str_c(CTL, SCZ, sep = "-"),
    pair_id = str_glue("Pair{str_pad(row_number(CTL), 2, pad = 0)}")
  )

hits <-
  c(
    "P38",
    "ERK",
    "JNK",
    "PIM"
  )

group_diffs <- samples_ids |>
  pmap(
    ~ generate_group_differences(..1, ..2,
                                 processed$dataset_max_exp,
                                 processed$scaled_models)
  ) |>
  set_names(samples_ids$comp_name)

plots <- group_diffs |>
  map( ~ generate_reverse_krsa_ggplot(.x, hits = hits)) |>
  map2(
    samples_ids$comp_name,
    ~ ggsave(
      str_glue("{.y}_reverse_krsa.png"),
      plot = .x,
      bg = "white",
      path = file.path("figures", "reverse_krsa"),
      width = 10,
      height = 8,
      units = "in",
      dpi = 300
    )
  )

hit_peptides <- KRSA_coverage_STK_PamChip_87102_v2 |>
  filter(Kin %in% hits) |>
  rename(Kinase = Kin,
         Peptide = Substrates)

combined_data <- group_diffs |>
  map(~ right_join(.x, hit_peptides, by = "Peptide")) |>
  map2(samples_ids$pair_id, ~ mutate(.x, Dataset = .y)) |>
  bind_rows() |>
  select(Dataset, Barcode, Peptide, Kinase, LFC, totalMeanLFC) |>
  write_csv(file.path("results", "filtered_reverse_krsa.csv"))


comparative_reverse_krsa <- function(kinase, dataset) {
  filtered <- dataset |>
    filter(Kinase == kinase)

  g <- ggplot(filtered, aes(x = Dataset, y = LFC))

  p <- g +
    geom_boxplot() +
    geom_jitter(width = 0.2, height = 0.2) +
    theme_minimal() +
    ggtitle(str_glue("Comparison of datasets for {kinase} family")) +
    scale_y_continuous(
      name = "Log_2 Fold Change",
      limits = c(-2.5, 2.5),
      breaks = seq(-2.5, 2.5, 0.5)
    ) +
    scale_x_discrete(name = "Dataset") +
    theme(plot.title = element_text(hjust = 0.5))


  p
}


hits |>
  map( ~ comparative_reverse_krsa(.x, combined_data)) |>
  walk(
    ~ ggsave(
      str_glue(
        "{unique(.x$data$Kinase)}_reverse_krsa_comparison.png"
      ),
      plot = .x,
      bg = "white",
      path = file.path("figures", "reverse_krsa"),
      width = 10,
      height = 8,
      units = "in",
      dpi = 300
    )
  )
