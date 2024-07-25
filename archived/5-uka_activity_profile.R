# Results from UKA
#

library(tidyverse)

mapping <- creedenzymatic::kinome_mp_file

erk_genes <-
  c("MEKK2",
    "MEKK3",
    "TPL2",
    "ERK5",
    "RAF1",
    "RAFA",
    "RAFB",
    "MEK1",
    "MEK2",
    "ERK1",
    "ERK2")
jnk_genes <-
  c(
    "MEKK1",
    "MEKK2",
    "MEKK3",
    "DLK",
    "MLK2",
    "TPL2",
    "ASK1",
    "TAK1",
    "TAO1",
    "TAO2",
    "MEK4",
    "MEK7",
    "JNK1",
    "JNK2",
    "JNK3"
  )
p38_genes <-
  c(
    "MEKK1",
    "MEKK2",
    "MEKK3",
    "DLK",
    "MLK2",
    "TPL2",
    "ASK1",
    "TAK1",
    "TAO1",
    "TAO2",
    "MEK3",
    "MEK6",
    "P38A",
    "P38B",
    "P38G",
    "P38D"
  )

genes <- read_csv("raw/Kinase-Expression.csv") |>
  drop_na(Name) |>
  pull(HGNC_Symbol)

uka_names <- mapping |>
  filter(hgnc_symbol %in% genes) |>
  pull(uka_id)

uka_results <- list.files("data/UKA_Reports/", full.names = TRUE) |>
  set_names(str_remove(list.files("data/UKA_Reports/"), "_.*")) |>
  map(~ read_tsv(.x)) |>
  map(~ select(.x, any_of(c("Kinase Name", "Mean Kinase Statistic")))) |>
  map2_dfr(str_remove(list.files("data/UKA_Reports/"), "_.*"), ~ mutate(.x, Dataset = .y)) |>
  filter(`Kinase Name` %in% uka_names) |>
  rename(Kinase = `Kinase Name`,
         Score = `Mean Kinase Statistic`)

g <- ggplot(uka_results, aes(x = Dataset, y = Kinase, fill = Score))

p <- g + geom_tile() +
  scale_fill_viridis_c()

