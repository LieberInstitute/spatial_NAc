library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library("tidyverse")

spe_in = here("processed-data", "05_harmony_BayesSpace", "spe_harmony.rds")
cluster_dir = here(
    "processed-data", "05_harmony_BayesSpace", "clusters_BayesSpace"
)
k2_cluster_name = "BayesSpace_harmony_k02"

#   Load SpatialExperiment and add k = 2 clusters
spe <- readRDS(spe_in)
cluster_import(
    spe,
    k2_cluster_name,
    cluster_dir = cluster_dir
)
