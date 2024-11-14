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

# Read in settings
spec <- matrix(
    c("cluster_dir", "d", 1, "character", "Specify the name of the cluster directory", 
      "cluster_file", "f", 1, "character", "Specify the name of the cluster file", 
       "agg_level", "a", 1, "character", "Aggregating at the donor or capture area level?", 
       "out_path", "o", 1, "character", "Output directory name?"),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)
#opt <- list()
#opt$cluster_dir <- "05_harmony_BayesSpace/05-BayesSpace_k_search/BayesSpace_harmony_k2"
#opt$cluster_file <- "clusters.csv"
#opt$agg_level <- "sample_id_original"
#opt$out_path <- "02_BayesSpace/pseudobulk_capture_area"

if(grepl("BayesSpace", opt$cluster_dir)){
    clust_dir <- unlist(strsplit(opt$cluster_dir, split = "/"))[3]
    clust_dir_n <- as.numeric(gsub("BayesSpace_harmony_k", "", clust_dir))
    clust_dir_n_nice <- sprintf("%02d", clust_dir_n)
    clust_dir_n <- paste0("BayesSpace_harmony_k", clust_dir_n_nice)
    clust_dir <- paste(unlist(strsplit(opt$cluster_dir, split = "/"))[1], unlist(strsplit(opt$cluster_dir, split = "/"))[2], 
    clust_dir_n, sep = "/")
    opt$cluster_dir <- clust_dir
}

print(opt$cluster_dir)
print(opt$cluster_file)
print(opt$agg_level)
print(opt$out_path)

# Read in the spatial data
raw_in_path <- here(
    "processed-data", "05_harmony_BayesSpace", "02-compute_QC_metrics", "spe_with_QC_metrics_hdf5"
)
spe <- loadHDF5SummarizedExperiment(raw_in_path)

# Specify the settings 
clusters_resDir <- here("processed-data", opt$cluster_dir)
clusters_resFile <- file.path(clusters_resDir, opt$cluster_file)
cluster_col <- gsub(".csv", "", opt$cluster_file)
cluster_col <- tolower(cluster_col)
if(grepl("BayesSpace", opt$cluster_dir)){
    spe_pseudo_path <- here("processed-data", "10_post_clustering_analysis", "01_pseudobulk_markers", opt$out_path, sprintf('spe_pseudo_%s.rds', unlist(strsplit(opt$cluster_dir, "/"))[3])) 
}
if(grepl("precast", opt$cluster_dir)){
    spe_pseudo_path <- here("processed-data", "10_post_clustering_analysis", "01_pseudobulk_markers", opt$out_path, sprintf('spe_pseudo_%s.rds', cluster_col)) 
}
# Add clustering results
#   Add a colData column for PRECAST results at this k value
if(grepl("BayesSpace", opt$cluster_dir)){
    spe$key <- paste0(spe$key, "_", spe$sample_id)
}

spe[[cluster_col]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(clusters_resFile), by = 'key') |>
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

#   Simplify colData to key, sample-level information
if(opt$agg_level == "sample_id_original"){
    colData(spe_pseudo) = colData(spe_pseudo)[
    , sort(c("sample_id", "sample_id_original", "slide_num", "in_tissue", "slide_num", "donor",  "Age", "Sex", cluster_col, "Diagnosis", "ncells"))
]
}
if(opt$agg_level == "sample_id"){
    colData(spe_pseudo) = colData(spe_pseudo)[
    , sort(c("sample_id", "slide_num", "in_tissue", "slide_num", "donor",  "Age", "Sex", cluster_col, "Diagnosis", "ncells"))
]
}

saveRDS(spe_pseudo, spe_pseudo_path)

session_info()

