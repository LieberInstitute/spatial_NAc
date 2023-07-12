

## Required libraries
library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library("BayesSpace")
library("Polychrome")

set.seed(20230712)

## Create output directories
dir_plots <- here::here("plots", "08_harmony_BayesSpace", opt$spetype)
dir_rdata <- here::here("processed-data", "08_harmony_BayesSpace", opt$spetype)
spe_in = file.path(dir_rdata, "spe_harmony.rds")

dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(dir_rdata, "clusters_BayesSpace"), showWarnings = FALSE)

## Load the data
spe <- readRDS(spe_in)

## Choose k
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))
k_nice <- sprintf("%02d", k)

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

message("Running spatialCluster()")
Sys.time()
spe <- spatialCluster(spe, use.dimred = "HARMONY", q = k)
Sys.time()

spe$bayesSpace_temp <- spe$spatial.cluster
bayesSpace_name <- paste0("BayesSpace_harmony_k", k_nice)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = file.path(dir_rdata, "clusters_BayesSpace")
)

## Visualize BayesSpace results
sample_ids <- unique(spe$sample_id)
cols <- Polychrome::palette36.colors(k)
names(cols) <- sort(unique(spe$spatial.cluster))

vis_grid_clus(
    spe = spe,
    clustervar = paste0("BayesSpace_harmony_k", k_nice),
    pdf_file = file.path(dir_plots, paste0("BayesSpace_harmony_k", k_nice, ".pdf")),
    sort_clust = FALSE,
    colors = cols,
    spatial = FALSE,
    point_size = 2
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
