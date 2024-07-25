# Reverse KRSA on the rats data

library(tidyverse)
library(KRSA)
library(gt)
library(knitr)

signal_file <- "raw/rats-haloperidol/rats-haloperidol_SigmBg.txt"
saturation_file <-
  "raw/rats-haloperidol/rats-haloperidol_Saturation.txt"

dataset_max <- krsa_read(signal_file, saturation_file) |>
  filter(SampleName == "Brain") |>
  krsa_qc_steps() |>
  mutate(
    Group = if_else(`Comment 1` == "CT", "CTL", "HPD"),
    SampleName = str_c(SampleName, Group, sep = "_")
  ) |>
  krsa_extractEndPointMaxExp(type = "STK")

dataset_all <- krsa_read(signal_file, saturation_file) |>
  filter(SampleName == "Brain") |>
  krsa_qc_steps() |>
  mutate(
    Group = if_else(`Comment 1` == "CT", "CTL", "HPD"),
    SampleName = str_c(SampleName, Group, sep = "_")
  ) |>
  krsa_extractEndPoint(type = "STK")

scaled_model <- dataset_max |>
  krsa_filter_lowPeps(threshold = 5) |>
  {
    \(x) krsa_scaleModel(dataset_all, x)
  }()


comparisons <- list(HPD = c("HPD", "CTL"))

selected_peptides_compared <-
  krsa_quick_filter(
    data = dataset_max,
    data2 = scaled_model$scaled,
    signal_threshold = 5,
    r2_threshold = 0.7,
    groups = comparisons$HPD
  )

hits <- c("JNK", "ERK", "P38")

hit_peptides <- KRSA_coverage_STK_PamChip_87102_v2 |>
  filter(Kin %in% hits) |>
  rename(Kinase = Kin,
         Peptide = Substrates)

group_diff_df <-
  krsa_group_diff(scaled_model$scaled,
                  comparisons$HPD,
                  selected_peptides_compared,
                  byChip = T)


combined <- group_diff_df |>
  right_join(hit_peptides, by = "Peptide")

directions <- combined |>
  mutate(Direction = case_when(
    LFC >= 0 ~ 1,
    LFC <= -0 ~ -1,
    TRUE ~ 0
  )) |>
  drop_na(LFC) |>
  group_by(Kinase) |>
  select(Direction) |>
  summarise(Kinase_Direction = mean(Direction)) |>
  mutate(Kinase_Direction = round(Kinase_Direction, 3))


fig <-
  krsa_reverse_krsa_plot(KRSA_coverage_STK_PamChip_87102_v2,
                         group_diff_df,
                         hits,
                         0.2,
                         F)

fig
