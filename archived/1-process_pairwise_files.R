# Process the raw paired data into a usable format

library(tidyverse)

annotation <- "annotation/sample_matching.csv" %>%
  read_csv(col_types = cols(.default = col_character())) %>%
  mutate(
    UKA_File_Path = file.path("raw", UKA_Dir, UKA_Subdir, UKA_File_Name),
    UKA_File_Path = ifelse(str_detect(UKA_File_Path, "NA"), NA, UKA_File_Path),
    KRSA_Full_File_Path = file.path("raw", KRSA_Dir, KRSA_Full_File_Name),
    KRSA_Full_File_Path = ifelse(
      str_detect(KRSA_Full_File_Path, "NA"),
      NA,
      KRSA_Full_File_Path
    ),
    KEA3_Peptide_File_Path = file.path("raw", KRSA_Dir, KRSA_Peptide_File_Name),
    KEA3_Peptide_File_Path = ifelse(
      str_detect(KEA3_Peptide_File_Path, "NA"),
      NA,
      KEA3_Peptide_File_Path
    )
  ) %>%
  pivot_longer(cols = contains("Path"),
               names_to = "Source",
               values_to = "SourcePath") %>%
  mutate(
    App = case_when(
      str_detect(Source, "KRSA") ~ "KRSA",
      str_detect(Source, "UKA") ~ "UKA",
      str_detect(Source, "KEA3") ~ "KEA3",
      TRUE ~ "NA"
    ),
    App = ifelse(App == "NA", NA, App),
    DestinationPath = file.path("data", str_glue("{App}_Reports"), str_glue("{Designation}_{str_remove(Source, '_Path')}.txt"))
  )

annotation %>%
  select(SourcePath, DestinationPath) %>%
  mutate(overwrite = TRUE) %>%
  rename(from = SourcePath,
         to = DestinationPath) %>%
  pwalk(file.copy)
