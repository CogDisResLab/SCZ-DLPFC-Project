# Run, store and analyze drugfindr data

suppressPackageStartupMessages({
  library(drugfindR)
  library(tidyverse)
})

target <- "MAPK11"

result <- investigateTarget(target, inputLib = "KD", outputLib = "CP", filterThreshold = 1)
