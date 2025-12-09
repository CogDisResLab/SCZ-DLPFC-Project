# Take the AssayQuant data and clean it up.

suppressPackageStartupMessages({
    library(tidyverse)
    library(readxl)
})

filepath <- file.path("raw", "Run 8_8.28.25_SCZ brain homogonate.xlsx")

well_layout <- read_excel(filepath,
    sheet = "Samples_Plate Map",
    range = "C1:O9",
    col_names = TRUE,
    col_types = "text"
) |>
    rename(Row = ...1) |>
    pivot_longer(-Row, names_to = "Column", values_to = "Sample") |>
    mutate(Coordinate = str_c(Row, Column), )

lysates <- read_excel(filepath, range = "B17:B20", col_names = "Sample", col_types = "text") |>
    mutate(Class = "LYS")

controls <- read_excel(filepath, range = "D17:E27", col_names = c("Sample", "Class"), col_types = "text")

scz <- read_excel(filepath, range = "F17:G36", col_names = c("Sample", "Class"), col_types = "text")

sample_mapping <- lysates |>
    bind_rows(controls) |>
    bind_rows(scz) |>
    group_by(Class) |>
    mutate(Class = str_c(Class, str_pad(seq_along(Class), width = 2, pad = "0")))

assay_control_mapping <- well_layout |>
    select(Sample) |>
    left_join(sample_mapping, by = "Sample") |>
    mutate(Class = case_when(
        str_detect(Sample, fixed("no lysate")) ~ "BLNK",
        str_detect(Sample, fixed("no reporter")) ~ "BLNK_LYS",
        str_detect(Sample, "5000ng.*ERK") ~ "LYS_INH",
        str_detect(Sample, "5000ng.*HI") ~ "LYS_HI",
        str_detect(Sample, "125\\s*nM") ~ "ERK_0.125nM",
        str_detect(Sample, "75\\s*nM") ~ "ERK_0.750nM",
        str_detect(Sample, "25\\s*nM") ~ "ERK_0.250nM",
        str_detect(Sample, "2nM.*INH") ~ "ERK_2.000nM_INH",
        str_detect(Sample, "2\\s*nM") ~ "ERK_2.000nM",
        str_detect(Sample, "1.5\\s*nM") ~ "ERK_1.500nM",
        str_detect(Sample, "1\\s*nM") ~ "ERK_1.000nM",
        str_detect(Sample, "5\\s*nM") ~ "ERK_0.500nM",
        is.na(Sample) ~ "EMPTY",
        .default = Class
    ))

combined <- well_layout |>
    inner_join(assay_control_mapping) |>
    write_csv(file.path("kinome_data", "erk_assayquant_plate_layout.csv"))
