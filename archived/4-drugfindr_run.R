# drugfindR Run for our chosen genes

# DrugFindR Run
# Running DrugFindR on JNK, ERK and P38

library(drugfindR)
library(tidyverse)
library(furrr)

kd_genes <- read_csv("ancillary/kd_filter.csv") |>
  pull(hgnc_symbol) |>
  unique()

oe_genes <- read_csv("ancillary/oe_filter.csv") |>
  pull(hgnc_symbol) |>
  unique()

kd_input_lib <- "KD"

oe_input_lib <- "OE"

output_lib <- "CP"

filter_threshold <- c(0.85, 1)

similarity_threshold <- c(0.2)

paired <- c(TRUE)

discordant <- c(TRUE, FALSE)

kd_matrix <- expand_grid(target = kd_genes, input_lib = kd_input_lib, output_lib,
                         filter_threshold, similarity_threshold, paired, discordant)

oe_matrix <- expand_grid(target = oe_genes, input_lib = oe_input_lib, output_lib,
                         filter_threshold, similarity_threshold, paired, discordant)

option_matrix <- bind_rows(kd_matrix, oe_matrix)

results <- option_matrix |>
  pmap(safely(investigate_target))

save.image(file="drugfindr.RData")
