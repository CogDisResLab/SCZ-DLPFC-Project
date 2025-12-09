# Load and clean the data from the assay quant results

suppressPackageStartupMessages({
    library(tidyverse)
    library(readxl)
    library(lubridate)
})

filepath <- file.path("raw", "Run 8_8.28.25_SCZ brain homogonate.xlsx")

dataset <- read_excel(filepath, sheet = "Raw Data") |>
    rename(Temperature = 2L) |>
    mutate(
        xConvertedTime = as_datetime(Time),
        xseconds = second(xConvertedTime),
        xminute = minute(xConvertedTime),
        xhour = hour(xConvertedTime),
        Elapsed = (xhour * 60L * 60L) + (xminute * 60L) + xseconds
    ) |>
    select(
        Elapsed, Temperature,
        starts_with("A"), starts_with("B"),
        starts_with("C"), starts_with("D"),
        starts_with("E"), starts_with("F"),
        starts_with("G"), starts_with("H")
    ) |>
    mutate(across(A1:H12, \(x) {
        str_x <- as.character(x)
        if_else(str_x == "OVRFLW", Inf, as.numeric(str_x))
    })) |>
    pivot_longer(c(-Temperature, -Elapsed), names_to = "Coordinate", values_to = "Reading")

map_filepath <- file.path("kinome_data", "erk_assayquant_plate_layout.csv")

plate_layout <- read_csv(map_filepath)

enriched_dataset <- dataset |>
    inner_join(plate_layout, by = "Coordinate", relationship = "many-to-many") |>
    select(
        ElapsedTime = Elapsed, Row, Column, Coordinate,
        Sample, Class, Temperature, Reading
    ) |>
    write_csv(file.path("kinome_data", "erk_assayquant_enriched_readings.csv"))
