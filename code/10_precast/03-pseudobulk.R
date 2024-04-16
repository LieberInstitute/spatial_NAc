library(here)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(getopt)

k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
spec <- matrix(
    c(
        "agg_level", "a", 1, "character", "Specify aggregation level"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

cluster_col = paste0('precast_k', k)

spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5"
)
precast_path = here(
    'processed-data', '10_precast', 'nnSVG_precast', sprintf('PRECAST_k%s.csv', k)
)

if(opt$agg_level == "sample_id"){
    spe_pseudo_path = here(
    'processed-data', '10_precast', 'nnSVG_precast', 'pseudobulk_donor', sprintf('spe_pseudo_%s.rds', cluster_col)
)
}else{
    if(opt$agg_level == "sample_id_original"){
        spe_pseudo_path = here(
    'processed-data', '10_precast', 'nnSVG_precast', 'pseudobulk_capture_area', sprintf('spe_pseudo_%s.rds', cluster_col)
)
    }else{
        stop("Invalid input")

    }
}

spe <- loadHDF5SummarizedExperiment(spe_dir)

#   Add a colData column for PRECAST results at this k value
spe[[cluster_col]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(precast_path), by = 'key') |>
    pull(cluster) |>
    as.factor()

# Based on the aggregation level perform pseudo-bulking
if(opt$agg_level == "sample_id"){
    spe_pseudo = registration_pseudobulk(
    spe,
    var_registration = cluster_col,
    var_sample_id = "sample_id",
    min_ncells = 10
)}

if(opt$agg_level == "sample_id_original"){
    spe_pseudo = registration_pseudobulk(
    spe,
    var_registration = cluster_col,
    var_sample_id = "sample_id_original",
    min_ncells = 10
)}



#   Simplify colData to key, sample-level information
colData(spe_pseudo) = colData(spe_pseudo)[
    , sort(c("Age", "Sex", "sample_id_original", cluster_col, "Diagnosis"))
]

saveRDS(spe_pseudo, spe_pseudo_path)

session_info()
