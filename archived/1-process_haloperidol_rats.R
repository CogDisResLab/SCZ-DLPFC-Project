# Process the haloperidol-treated rats file

library(tidyverse)

source_prefix <- file.path("raw", "rats-haloperidol")
dest_prefix <- file.path("data", "Rats-Data")


files <- list.files(source_prefix, "HALOPERIDOL")

annotation <- files |>
  str_match("(HALOPERIDOL)_(\\w{3,4}).txt") |>
  as.data.frame() |>
  set_names(c("file", "treatment", "tool")) |>
  mutate(out_file = str_glue("Rats_HPD_{tool}.csv"),
         data = map(file.path(source_prefix, file), ~ read_tsv(.x)),
         out = map2(data, file.path(dest_prefix, out_file), ~ write_csv(.x, .y)))
