# Generate guided network for selected pairs

library(tidyverse)
library(KINNET)
library(bnlearn)

selected_samples <-
  c("R1C2P3",
    "R2C2P4",
    "R3C3P5",
    "R4C1P2",
    "R1C1P2",
    "R2C2P3",
    "R3C1P1",
    "R3C2P3")

stk_data <- KINNET::PamchipData_STK("raw/KINNET-Input/Run04_KRSAsOutput_Median_SigmBg_230705171954.txt")

filtered_peptides <- stk_data |>
  filter_peptides(threshold = 0.5)

ctl_subset <- subset_data(stk_data, filtered_peptides$CTL_882, "CTL_882")
scz_subset <- subset_data(stk_data, filtered_peptides$SCZ_585, "SCZ_585")

ctl_model <- make_model(ctl_subset, 1000, threshold = 0)
scz_model <- make_model(scz_subset, 1000, threshold = 0)

saveRDS(list(filtered_peptides, ctl_subset, ctl_model, scz_subset, scz_model), file = "depot/R4C1P2-Model-Data.RDS")

ctl_averaged <- averaged.network(ctl_model$strength_net, 0.7)
scz_averaged <- averaged.network(scz_model$strength_net, 0.7)

ctl_kinased <- assign_kinases_guided(ctl_averaged, "STK", guided = c("SGK1"))
scz_kinased <- assign_kinases_guided(scz_averaged, "STK", guided = c("SGK1"))

pdf("figures/R4C1P2-Comparison.pdf")
compare_kinased_graphs(ctl_kinased, scz_kinased, "CTL_832", "SCZ_585", render = TRUE)
dev.off()
