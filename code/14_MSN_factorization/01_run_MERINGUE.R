library(MERINGUE)
library(sessioninfo)
library(here)
library(HDF5Array)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)

# Set directories
plotDir <- here("plots", "14_MSN_factorization", "01_meringue")
resDir <- here("processed-data", "14_MSN_factorization", "01_meringue")
# Read in the processed SPE data
spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_dimRed_hdf5"
)
spe <- loadHDF5SummarizedExperiment(spe_dir)
clusters_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters", "precast_clusters.csv")
spe[["spatial_domains"]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(clusters_path), by = 'key') |>
    pull(cluster) |>
    as.factor()

# Remove spots with no assigned clusters
spe <- spe[ ,!is.na(spe[["spatial_domains"]])]
#pdf(file.path(plotDir, "preliminary_clustering_visualization.pdf"), height = 10, width = 10)
#vis_clus(spe, sampleid = "Br6522", clustervar = "spatial_domains", spatial = TRUE, is_stitched = TRUE) +
            #   Increase size of colored dots in legend
#            guides(fill = guide_legend(override.aes = list(size = 5)))
#vis_clus(spe, sampleid = "Br6432", clustervar = "spatial_domains", spatial = TRUE, is_stitched = TRUE) +
            #   Increase size of colored dots in legend
#            guides(fill = guide_legend(override.aes = list(size = 5)))
#vis_clus(spe, sampleid = "Br3942", clustervar = "spatial_domains", spatial = TRUE, is_stitched = TRUE) +
            #   Increase size of colored dots in legend
#            guides(fill = guide_legend(override.aes = list(size = 5)))
#dev.off()

# Format the spatialCoords to work with the gene expression data
rownames(spatialCoords(spe)) = colnames(spe)

# Subset to a particular sample
sample_id <- levels(spe$sample_id)[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]
spe <- spe[ ,spe$sample_id == sample_id]
spe <- spe[ ,grepl("MSN", spe[["spatial_domains"]])]

#   Grab a corner (1% of the total area) for testing
#a = spatialCoords(spe)[,1] > 0.1 * min(spatialCoords(spe)[, 1]) + 0.8 * max(spatialCoords(spe)[, 1])
#b = spatialCoords(spe)[,2] < 0.1 * min(spatialCoords(spe)[, 2]) + 0.9 * max(spatialCoords(spe)[, 2])
#spe = spe[, a& b]

#pdf(file.path(plotDir, "test_spots.pdf"), height = 10, width = 10)
#vis_clus(spe, sampleid = sample_id, clustervar = "spatial_domains", spatial = TRUE, is_stitched = TRUE) +
            #   Increase size of colored dots in legend
#            guides(fill = guide_legend(override.aes = list(size = 5)))
#dev.off()

# Run MERINGUE
# Get spatial neighbors for the positions
message(Sys.time(), ' | Forming network of spatial neighbors')
W <- getSpatialNeighbors(spatialCoords(spe), filterDist = 500)

usePreComputedI <- TRUE
if(usePreComputedI){
    I <- readRDS(file.path(resDir, paste0("I_", sample_id, ".rds")))
    mat <- logcounts(spe)
    rownames(mat) <- make.names(rownames(mat), unique = TRUE)
}else{
    png(
    file.path(
        plotDir, sprintf('neighbor_network_%s.png', sample_id)
    ),
    width = 1500, height = 1500
    )
    plotNetwork(spatialCoords(spe), W)
    dev.off()

    # Identify sigificantly spatially auto-correlated genes
    start_time <- Sys.time()
    mat <- logcounts(spe)
    rownames(mat) <- make.names(rownames(mat), unique = TRUE)
    I <- getSpatialPatterns(mat, W)
    end_time <- Sys.time()
    print(end_time - start_time)

    # Save raw results
    saveRDS(I, file.path(resDir, paste0("I_", sample_id, ".rds")))
}

results.filter <- filterSpatialPatterns(mat = mat,
                                        I = I,
                                        w = W,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = 0.05,
                                        verbose = TRUE)
saveRDS(results.filter, file.path(resDir, paste0("results_filtered_", sample_id, ".rds")))
# Compute spatial cross correlation matrix
scc <- spatialCrossCorMatrix(mat = as.matrix(mat[results.filter,]), 
                             weight = W)
saveRDS(scc, file.path(resDir, paste0("scc_", sample_id, ".rds")))