# Significant Peptide Selection by PLSDA

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
  purrr::imap_dfr( ~ mutate(.x, Run = .y)) |>
  filter(str_detect(Peptide, "^#REF.*", negate = TRUE))



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

batch_corrected <- batched |>
  t()

grouping <- heatmap_annotation$Group |> as.factor()

## Select Parameters for optimal PLS-DA

### Select n of components

comparative_plsda <- splsda(batch_corrected, grouping, ncomp = 15)

perf_comparative_plsda <- perf(comparative_plsda, validation = "Mfold", folds = 5,
                              progressBar = TRUE, nrepeat = 50)

list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX

tune_plsda <- tune.splsda(batch_corrected, grouping, ncomp = 5, validation = 'Mfold',
                          folds = 5, dist = 'max.dist', progressBar = TRUE,
                          measure = "BER", test.keepX = list.keepX,
                          nrepeat = 50)

error <- tune_plsda$error.rate
ncomp <- tune_plsda$choice.ncomp$ncomp

select_keepX <- tune_plsda$choice.keepX[1:ncomp]


plsda_result <- splsda(batch_corrected, grouping, ncomp = ncomp, keepX = select_keepX)

background <- background.predict(plsda_result, comp.predicted = ncomp)

plotIndiv(plsda_result, comp = 1:ncomp, group = grouping, ind.names = FALSE, legend = TRUE, background = background, ellipse = TRUE)


selected_peptides_1 <- rownames(selectVar(plsda_result, comp = 1)$value)
selected_peptides_2 <- rownames(selectVar(plsda_result, comp = 2)$value)

all_peptides <- union(selected_peptides_1, selected_peptides_2)
