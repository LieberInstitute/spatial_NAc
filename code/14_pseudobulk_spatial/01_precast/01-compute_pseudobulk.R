library(here)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(getopt)
library(edgeR)
library(scater)
library(scran)
library(dplyr)
library(gridExtra)
library(ggforce)
library(pheatmap)


k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
spec <- matrix(
    c("agg_level", "a", 1, "character", "Specify whether to aggregate at the donor or capture area level?"),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

cluster_col = paste0('precast_k', k)

raw_in_path <- here(
    "processed-data", "05_harmony_BayesSpace", "02-compute_QC_metrics", "spe_with_QC_metrics_hdf5"
)
precast_path = here(
    'processed-data', '10_precast', 'nnSVG_precast', sprintf('PRECAST_k%s.csv', k)
)

if(opt$agg_level == "sample_id_original"){
    spe_pseudo_path = here(
    'processed-data', '14_pseudobulk_spatial', '01_precast', 'pseudobulk_capture_area', sprintf('spe_pseudo_%s.rds', cluster_col)
)
}else{
    if(opt$agg_level == "sample_id"){
        spe_pseudo_path = here(
    'processed-data', '14_pseudobulk_spatial', '01_precast', 'pseudobulk_donor', sprintf('spe_pseudo_%s.rds', cluster_col)
)
    }else{
        stop("Incorrect input for aggregation level")
    }
}

spe <- loadHDF5SummarizedExperiment(raw_in_path)

#   Add a colData column for PRECAST results at this k value
spe[[cluster_col]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(precast_path), by = 'key') |>
    pull(cluster) |>
    as.factor()

# Remove spots with no PRECAST output
spe <- spe[ ,!is.na(spe[[cluster_col]])]

if(opt$agg_level == "sample_id_original"){
    spe_pseudo <- aggregateAcrossCells(
    spe, DataFrame(cluster = spe[[cluster_col]], sample_id = colData(spe)$sample_id_original))
    colnames(spe_pseudo) <- paste0(spe_pseudo$sample_id_original, "_", spe_pseudo[[cluster_col]])
}else{
    spe_pseudo <- aggregateAcrossCells(
    spe, DataFrame(cluster = spe[[cluster_col]], sample_id = colData(spe)$sample_id))
    colnames(spe_pseudo) <- paste0(spe_pseudo$sample_id, "_", spe_pseudo[[cluster_col]])
}
spe_pseudo = registration_pseudobulk(
spe,
var_registration = cluster_col,
var_sample_id = "sample_id_original",
min_ncells = 0)

    dim(spe_pseudo)
colnames(spe_pseudo) <-
            paste0(
                spe_pseudo$registration_sample_id,
                "_",
                spe_pseudo$registration_variable
            )
#   Simplify colData to key, sample-level information
colData(spe_pseudo) = colData(spe_pseudo)[
    , sort(c("sample_id", "sample_id_original", "slide_num", "in_tissue", "slide_num", "donor",  "Age", "Sex", cluster_col, "Diagnosis", "ncells"))
]

saveRDS(spe_pseudo, spe_pseudo_path)

session_info()
