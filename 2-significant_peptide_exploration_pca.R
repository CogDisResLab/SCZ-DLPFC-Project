# Significant Peptide Selection by PCA/PLSDA/Threshold

library(mixOmics)
library(KRSA)
library(sva)
library(tidyverse)

qc_peptides <- function(run,
                        chip_type = "STK",
                        index = NULL) {
  data <- krsa_read(run$signal, run$saturation) |>
    krsa_qc_steps() |>
    mutate(Group = str_extract(SampleName, "^\\w{3}"),
           Run = index) |>
    filter(Group %in% c("SCZ", "CTL"))

  data_end_max <- krsa_extractEndPointMaxExp(data, "STK")
  data_end <- krsa_extractEndPoint(data, "STK")

  filtered_low_peps <- krsa_filter_lowPeps(data_end_max, 5)

  scaled_model <- krsa_scaleModel(data_end, filtered_low_peps)

  model <- scaled_model

  model
}


run1_files <- list(signal = "raw/Paired-Raw-Data/run1/_KRSAsfiles_Median_SigmBg_210518210643.txt",
                   saturation = "raw/Paired-Raw-Data/run1/_KRSAsfiles_Signal_Saturation_210518210643.txt")

run2_files <- list(signal = "raw/Paired-Raw-Data/run2/_KRSA_files_Median_SigmBg_210520201325.txt",
                   saturation = "raw/Paired-Raw-Data/run2/_KRSA_files_Signal_Saturation_210520201326.txt")

run3_files <- list(signal = "raw/Paired-Raw-Data/run3/_KRSA_files_Median_SigmBg_210519103637.txt",
                   saturation = "raw/Paired-Raw-Data/run3/_KRSA_files_Signal_Saturation_210519103637.txt")

run4_files <- list(signal = "raw/Paired-Raw-Data/run4/_KRSA_files_Median_SigmBg_210605104306.txt",
                   saturation = "raw/Paired-Raw-Data/run4/_KRSA_files_Signal_Saturation_210605104307.txt")

output <- list(run1_files, run2_files, run3_files, run4_files) |>
  imap( ~ qc_peptides(.x, index = .y))


all_data <- output |>
  purrr::map( ~ pluck(.x, "normalized")) |>
  purrr::imap_dfr( ~ mutate(.x, Run = .y))


heatmap_data <- all_data |>
  select(SampleName, Peptide, slope) |>
  pivot_wider(names_from = SampleName,
              values_from = slope,
              values_fill = 0) |>
  column_to_rownames("Peptide") |>
  as.matrix()

heatmap_annotation <- all_data |>
  select(SampleName, Group, Run) |>
  distinct() |>
  column_to_rownames("SampleName")

batched <- ComBat(heatmap_data, batch = heatmap_annotation$Run)

pca <- prcomp(t(heatmap_data), scale = TRUE)

pca.data <- data.frame(
  Sample = rownames(pca$x),
  X = pca$x[, 1],
  Y = pca$x[, 2],
  Run = factor(heatmap_annotation[rownames(pca$x), ]$Run),
  Group = heatmap_annotation[rownames(pca$x), ]$Group
)
pca.data

g <-
  ggplot(data = pca.data, aes(
    x = X,
    y = Y,
    color = Run,
    shape = Group
  ))

p <- g + geom_point() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("Batch-Corrected PCA of DLPFC SCZ Samples") +
  theme(legend.position = "bottom")

ggsave("figures/PCA-Classification-Results.png", plot = p, width = 10, height = 10, units = "in")
