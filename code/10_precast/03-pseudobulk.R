library(here)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)
library(SpatialExperiment)
library(getopt)
library(spatialLIBD)

spec <- matrix(
    c("k", "k", 1, "integer", "Number of PRECAST clusters"),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec = spec)

spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "spe_filtered_hdf5"
)
precast_path = here(
    'processed-data', '10_precast', sprintf('PRECAST_k%s.csv', opt$k)
)

spe <- loadHDF5SummarizedExperiment(spe_dir)

#   Add a colData column for PRECAST results at this k value
cluster_col = paste0('precast_k', opt$k)
spe[[cluster_col]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(precast_path), by = 'key') |>
    pull(cluster) |>
    as.factor()

spe_pseudo = registration_pseudobulk(
    spe,
    var_registration = cluster_col,
    var_sample_id = "sample_id",
    min_ncells = 10
)

#   Simplify colData to key, sample-level information
colData(spe_pseudo) = colData(spe_pseudo)[
    , sort(c("Age", "Sex", "sample_id", cluster_col, "Diagnosis"))
]

saveRDS(
    spe_pseudo,
    file.path(dirname(precast_path), sprintf('spe_pseudo_%s.rds', cluster_col))
)

session_info()
