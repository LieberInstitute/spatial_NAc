library(spacexr)
library(Matrix)
library(SingleCellExperiment)
library(here)
library(scran)
library(scater)
library(SpatialExperiment)
library(spatialLIBD)
library(spatialNAcUtils)
library(HDF5Array)
library(ggplot2)
library(patchwork)
library(getopt)

# Read the myRCTD object
spec <- matrix(
    c(
        "marker_genes", "m", 1, "logical", "Use only marker genes from single cell labels?"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)
spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5")
spe <- loadHDF5SummarizedExperiment(spe_dir)
sample_id <- levels(spe$sample_id)[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]

RCTD_dir <- here::here("processed-data", "08_spot_deconvo", "01_RCTD", sample_id)
if(opt$marker_genes){
    myRCTD = readRDS(here(RCTD_dir,"myRCTD_markers.rds"))
}else{
    myRCTD = readRDS(here(RCTD_dir,"myRCTD.rds"))
}

celltypes = levels(myRCTD@reference@cell_types)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'multi')
if(opt$marker_genes){
    saveRDS(myRCTD, file.path(RCTD_dir, "results_RCTD_markers.rds"))
}else{
    saveRDS(myRCTD, file.path(RCTD_dir, "results_RCTD.rds"))
}