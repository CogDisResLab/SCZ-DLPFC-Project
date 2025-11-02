# Extract and process the metadata for the SCZ subjects

suppressPackageStartupMessages({
    library(readxl)
    library(tidyverse)
})

metadata_file <- file.path("raw", "ERK-AQ-AD and SCZ Brains List.xlsx")

coords_table <- list(
    table1 = list(left_corner = "A14", right_corner = "J21"),
    table2 = list(left_corner = "A26", right_corner = "J37")
)

metadata_table1 <- read_excel(metadata_file, range = "A14:J21", na = c("", "N")) |>
    filter(!is.na(Age)) |>
    rename(Subject = 1L, Slide_Num = 2L) |>
    select(Subject, Slide_Num, Age, Race, pH, PMI, RIN, From) |>
    mutate(Gender = "M")

metadata_table2 <- read_excel(metadata_file, range = "A26:J37", na = c("", "N", "n/a")) |>
    filter(!is.na(Age)) |>
    rename(Subject = 1L, Slide_Num = 2L, PMI = 6, Laterality = 9) |>
    select(Subject, Slide_Num, Age, Race, pH, PMI, RIN, From) |>
    mutate(Gender = "F")

metadata <- bind_rows(metadata_table1, metadata_table2) |>
    write_csv(file.path("kinome_data", "erk_assayquant_subjects.csv"))
