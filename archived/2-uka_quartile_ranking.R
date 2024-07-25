# Create Quartile Rankings from UKA Data

library(tidyverse)

gendered <- read_csv("annotation/sample_matching.csv", col_types = cols(.default = col_character())) %>%
  select(Designation, Gender) %>%
  group_by(Gender) %>%
  group_split() %>%
  set_names(c("F", "M")) %>%
  map(~ ungroup(.x)) %>%
  map(~ pull(.x, Designation))

files <- list.files("data/UKA_Reports/", full.names = TRUE)

filenames <- files %>%
  map_chr(basename) %>%
  map_chr( ~ str_extract(.x, "R\\dC\\dP\\d"))

col_spec <- cols(
  .default = col_skip(),
  `Kinase Name` = col_character(),
  `Median Final score` = col_double())

datasets <- files %>%
  map(~ read_tsv(.x, col_types = col_spec)) %>%
  set_names(filenames) %>%
  map(~ rename(.x, Score = `Median Final score`, Kinase = `Kinase Name`)) %>%
  map(~ arrange(.x, desc(abs(Score)))) %>%
  map2_dfr(filenames,
           ~ mutate(
             .x,
             dataset = .y,
             rank = row_number(desc(Score)),
             quartile = ntile(desc(Score), 4)
           )) %>%
  write_csv("results/UKA_quartile_rankings_all.csv") %>%
  select(-Score, -rank) %>%
  pivot_wider(names_from = dataset,
              values_from = quartile,
              values_fill = 100) %>%
  write_csv("results/UKA_quartile_rankings_matrix_all.csv") %>%
  pivot_longer(cols = where(is.numeric),
               names_to = "Dataset",
               values_to = "QuartileRank") %>%
  mutate(Quartile = factor(str_c("Q", QuartileRank)),
         Gender = case_when(
           Dataset %in% gendered$F ~ "F",
           Dataset %in% gendered$M ~ "M"
         ))

sum_ordered_kinases <- datasets %>%
  group_by(Kinase) %>%
  summarise(Total_Rank = mean(QuartileRank)) %>%
  ungroup() %>%
  arrange(Total_Rank) %>%
  write_csv("results/UKA_quartile_rankings_aggregated_ordered_all.csv")


sum_ordered_kinases_male <- datasets %>%
  filter(Dataset %in% gendered$M) %>%
  group_by(Kinase) %>%
  summarise(Total_Rank = mean(QuartileRank)) %>%
  ungroup() %>%
  arrange(Total_Rank) %>%
  write_csv("results/UKA_quartile_rankings_aggregated_ordered_male.csv")

sum_ordered_kinases_female <- datasets %>%
  filter(Dataset %in% gendered$F) %>%
  group_by(Kinase) %>%
  summarise(Total_Rank = mean(QuartileRank)) %>%
  ungroup() %>%
  arrange(Total_Rank) %>%
  write_csv("results/UKA_quartile_rankings_aggregated_ordered_female.csv")
