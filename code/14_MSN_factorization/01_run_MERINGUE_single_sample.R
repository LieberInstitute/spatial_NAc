library(MERINGUE)
library(sessioninfo)
library(here)
library(HDF5Array)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)
library(Rcpp)

# Compile the C++ function
sourceCpp(here("code", "14_MSN_factorization", "spatial_fast.cpp"))

# Define your wrapper
spatialCrossCorMatrix_fast <- function(mat, weight) {
 # Expect: mat = genes x cells
  #         weight = cells x cells

  if (nrow(weight) != ncol(weight)) {
    stop("'weight' must be a square matrix (cells × cells)")
  }

  if (ncol(mat) != nrow(weight)) {
    stop("'mat' must be genes × cells, where ncol(mat) == nrow(weight)")
  }

  if (!all(rownames(weight) %in% colnames(mat))) {
    stop("Column names of 'mat' (cells) must match rownames of 'weight'")
  }

  # Reorder columns of mat to match weight
  mat <- mat[, rownames(weight)]

  # Convert weight to sparse if necessary
  if (!inherits(weight, "dgCMatrix")) {
    weight <- as(weight, "dgCMatrix")
  }

  # Compute spatial cross-correlation
  scc <- spatialCrossCorMatrix_sparse_C(as.matrix(mat), weight)

  # Name output
  rownames(scc) <- rownames(mat)
  colnames(scc) <- rownames(mat)

  return(scc)
}

# Set directories
plotDir <- here("plots", "14_MSN_factorization", "01_meringue")
resDir <- here("processed-data", "14_MSN_factorization", "01_meringue")
# Read in the processed SPE data
spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_dimRed_hdf5"
)
spe <- loadHDF5SummarizedExperiment(spe_dir)


# Add clustering analysis from PRECAST
clusters_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters", "precast_clusters.csv")
spe[["spatial_domains"]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(clusters_path), by = 'key') |>
    pull(cluster) |>
    as.factor()

# Remove spots with no assigned clusters
spe <- spe[ ,!is.na(spe[["spatial_domains"]])]

# Format the spatialCoords to work with the gene expression data
rownames(spatialCoords(spe)) = colnames(spe)

# Subset to a particular sample if not transforming the row and column
sample_id <- levels(spe$sample_id)[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]
spe <- spe[ ,spe$sample_id == sample_id]
spe <- spe[ ,grepl("MSN", spe[["spatial_domains"]])]

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

usePreComputedresults <- TRUE
if(usePreComputedresults){
    results.filter <- readRDS(file.path(resDir, paste0("results_filtered_", sample_id, ".rds")))
}else{
    results.filter <- filterSpatialPatterns(mat = mat,
                                        I = I,
                                        w = W,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = 0.05,
                                        verbose = TRUE)
    saveRDS(results.filter, file.path(resDir, paste0("results_filtered_", sample_id, ".rds")))
}


# Compute spatial cross correlation matrix
scc <- spatialCrossCorMatrix_fast(mat = as.matrix(mat[results.filter,]), 
                             weight = W)
saveRDS(scc, file.path(resDir, paste0("scc_", sample_id, ".rds")))