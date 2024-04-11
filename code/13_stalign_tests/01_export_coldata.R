library(here)
library(tidyverse)
library(SpatialExperiment)

spe_path = here(
    'processed-data', '05_harmony_BayesSpace', '01-build_spe', 'spe_raw.rds'
)
out_path = here(
    'processed-data', '13_stalign_tests', 'raw_coldata.csv'
)

spe = readRDS(spe_path)

colData(spe) |>
    as_tibble() |>
    write_csv(out_path)
