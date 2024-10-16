## Required libraries
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("BayesSpace")
    library("RColorBrewer")
    library("ggplot2")
    library("gridExtra")
    library("Polychrome")
})

set.seed(20230712)

## Choose k
k <- as.numeric(
    #   Only one of these environment variables will be defined, so grab the
    #   defined one (handle SGE or SLURM)
    paste0(Sys.getenv("SLURM_ARRAY_TASK_ID"), Sys.getenv("SGE_TASK_ID"))
)
k_nice <- sprintf("%02d", k)

## Create output directories
dir_plots <- here("plots", "05_harmony_BayesSpace", "05-BayesSpace_k_search", k_nice)
dir_rdata <- here("processed-data", "05_harmony_BayesSpace", "05-BayesSpace_k_search")
spe_in = file.path(here("processed-data", "05_harmony_BayesSpace", "04-preprocess_and_harmony") ,"spe_harmony.rds")

dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(dir_rdata, "clusters_BayesSpace"), showWarnings = FALSE)

## Load the data
spe <- readRDS(spe_in)

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)
colData(spe)$row <- spe$array_row
colData(spe)$col <- spe$array_col

message("Running spatialCluster()")
Sys.time()
spe <- spatialCluster(spe, use.dimred = "HARMONY", q = k, nrep = 10000)
Sys.time()

spe$bayesSpace_temp <- as.factor(spe$spatial.cluster)
bayesSpace_name <- paste0("BayesSpace_harmony_k", k_nice)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = file.path(dir_rdata, bayesSpace_name)
)

## Visualize BayesSpace results
#sample_ids <- unique(spe$sample_id)
#cols <- Polychrome::palette36.colors(k)
#names(cols) <- sort(unique(spe[[bayesSpace_name]]))

#   Use 'vis_grid_clus' to preserve all spots (including overlaps)
#p_list = vis_grid_clus(
#    spe = spe,
#    clustervar = bayesSpace_name,
#    sort_clust = FALSE,
#    colors = cols,
#    spatial = FALSE,
#    point_size = 1,
#    auto_crop = FALSE,
#    return_plots = TRUE
#)
#pdf(file.path(dir_plots, paste0(bayesSpace_name, "_raw.pdf")))
#print(p_list)
#dev.off()

#   Use 'spot_plot', which takes one spot in case of overlaps
#p_list = list()
#for (sample_id in sample_ids) {
#    p_list[[sample_id]] = spot_plot(
#        spe, sample_id = sample_id, var_name = bayesSpace_name, colors = cols
#    )
#}
#pdf(file.path(dir_plots, paste0(bayesSpace_name, "_fit.pdf")))
#print(p_list)
#dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
