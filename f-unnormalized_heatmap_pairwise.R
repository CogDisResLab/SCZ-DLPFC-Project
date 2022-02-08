# Giant Heatmap of all peptides across all pairs - Without Normalization

library(tidyverse)
library(KRSA)

r1 <- krsa_read("raw/Paired-Raw-Data/run1/_KRSAsfiles_Median_SigmBg_210518210643.txt",
                "raw/Paired-Raw-Data/run1/_KRSAsfiles_Signal_Saturation_210518210643.txt")

r2 <- krsa_read("raw/Paired-Raw-Data/run2/_KRSA_files_Median_SigmBg_210520201325.txt",
                "raw/Paired-Raw-Data/run2/_KRSA_files_Signal_Saturation_210520201326.txt")

r3 <- krsa_read("raw/Paired-Raw-Data/run3/_KRSA_files_Median_SigmBg_210519103637.txt",
                "raw/Paired-Raw-Data/run3/_KRSA_files_Signal_Saturation_210519103637.txt")

r4 <- krsa_read("raw/Paired-Raw-Data/run4/_KRSA_files_Median_SigmBg_210605104306.txt",
                "raw/Paired-Raw-Data/run4/_KRSA_files_Signal_Saturation_210605104307.txt") |>
  filter(str_detect(SampleName, "SCZ") | str_detect(SampleName, "CTL"))


r <- bind_rows(r1, r2, r3, r4) |>
  krsa_qc_steps() |>
  mutate(Group = str_extract(SampleName, "^\\w{3}"))

r_max_exp <- krsa_extractEndPointMaxExp(r, "STK")

r_max <- krsa_extractEndPoint(r, "STK")

qc_passed_peps <- krsa_filter_lowPeps(r_max_exp, 5)

r_scaled <- krsa_scaleModel(r_max, qc_passed_peps)

final_passed_peps <- krsa_filter_nonLinear(r_scaled$scaled, 0.8) |>
  krsa_filter_ref_pep()

png(file.path("figures", "unnormalized_unscaled-all_pairs_heatmap.png"), width = 11, height = 8.5, units = "in", res = 300)
krsa_heatmap(r_scaled$scaled, final_passed_peps, scale = "none")
dev.off()

png(file.path("figures", "unnormalized_scaled-all_pairs_heatmap.png"), width = 11, height = 8.5, units = "in", res = 300)
krsa_heatmap(r_scaled$scaled, final_passed_peps, scale = "row")
dev.off()
