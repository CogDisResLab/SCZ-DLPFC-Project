# Clustering samples

library(tidyverse)
library(sva)
library(KRSA)
library(pheatmap)

qc_peptides <- function(run, chip_type = "STK", index = NULL) {

  data <- krsa_read(run$signal, run$saturation) |>
    krsa_qc_steps() |>
    mutate(Group = str_extract(SampleName, "^\\w{3}"),
           Run = index) |>
    filter(Group %in% c("SCZ", "CTL"))

  data_end_max <- krsa_extractEndPointMaxExp(data, "STK")
  data_end <- krsa_extractEndPoint(data, "STK")

  filtered_low_peps <- krsa_filter_lowPeps(data_end_max, 5)

  scaled_model <- krsa_scaleModel(data_end, filtered_low_peps)

  filtered_final <- scaled_model |>
    pluck("scaled") |>
    krsa_filter_nonLinear(0.8) |>
    krsa_filter_ref_pep()

  out <- list(data = data,
              peptides = filtered_final,
              model = scaled_model)

  out
}


run1_files <- list(
  signal = "raw/Paired-Raw-Data/run1/_KRSAsfiles_Median_SigmBg_210518210643.txt",
  saturation = "raw/Paired-Raw-Data/run1/_KRSAsfiles_Signal_Saturation_210518210643.txt"
)

run2_files <- list(
  signal = "raw/Paired-Raw-Data/run2/_KRSA_files_Median_SigmBg_210520201325.txt",
  saturation = "raw/Paired-Raw-Data/run2/_KRSA_files_Signal_Saturation_210520201326.txt"
)

run3_files <- list(
  signal = "raw/Paired-Raw-Data/run3/_KRSA_files_Median_SigmBg_210519103637.txt",
  saturation = "raw/Paired-Raw-Data/run3/_KRSA_files_Signal_Saturation_210519103637.txt"
)

run4_files <- list(
  signal = "raw/Paired-Raw-Data/run4/_KRSA_files_Median_SigmBg_210605104306.txt",
  saturation = "raw/Paired-Raw-Data/run4/_KRSA_files_Signal_Saturation_210605104307.txt"
)

output <- list(run1_files, run2_files, run3_files, run4_files) |>
  imap(~ qc_peptides(.x, index = .y))


selected_peptides <- output |>
  map(~ pluck(.x, "peptides")) |>
  reduce(union)

all_data <- output |>
  map(~ pluck(.x, "model")) |>
  map(~ pluck(.x, "normalized")) |>
  imap_dfr(~ mutate(.x, Run = .y))


heatmap_data <- all_data |>
  filter(Peptide %in% selected_peptides) |>
  select(SampleName, Peptide, slope) |>
  pivot_wider(names_from = SampleName, values_from = slope, values_fill = 0) |>
  column_to_rownames("Peptide") |>
  as.matrix()

heatmap_annotation <- all_data |>
  select(SampleName, Group, Run) |>
  distinct() |>
  column_to_rownames("SampleName")

batched <- ComBat(heatmap_data, batch = heatmap_annotation$Run)

distance <- dist(scale(t(batched)))

p <- pheatmap(batched, clustering_distance_cols = distance, annotation_col = heatmap_annotation, scale = "row", color = colorRampPalette(c("yellow", "white", "red"))(n = 50))

ggsave("figures/all_sample_clustering.png", bg = "white", width = 7, height = 10, units = "in")
