library(tidyverse)
library(here)
library(SpatialExperiment)
library(SingleCellExperiment)
library(zellkonverter)
library(HDF5Array)

spe_path <- here("processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5")
output_dir <- here("processed-data","22_gene_risk_LR_analysis", "03_tensor_cell2cell", "input_files")

# Load objects
spe <- loadHDF5SummarizedExperiment(spe_path)

# Add the metadata
opt <- list()
opt$marker_genes <- TRUE
spe_Br2743 <- mirrorObject(spe, sample_id = "Br2743", image_id = "lowres", axis = "v")
spe_Br8492 <- mirrorObject(spe, sample_id = "Br8492", image_id = "lowres", axis = "v")
spe_Br8325 <- mirrorObject(spe, sample_id = "Br8325", image_id = "lowres", axis = "v")
spe_Br3942 <- mirrorObject(spe, sample_id = "Br3942", image_id = "lowres", axis = "v")

spe <- spe[ ,!spe$sample_id %in% c("Br2743", "Br8492", "Br8325", "Br3942")]
spe <- cbind(spe, spe_Br2743)
spe <- cbind(spe, spe_Br8492)
spe <- cbind(spe, spe_Br8325)
spe <- cbind(spe, spe_Br3942)

sample_ids <- levels(spe$sample_id)
RCTD_list <- lapply(sample_ids, function(isample){
    RCTD_dir <- here::here("processed-data", "08_spot_deconvo", "01_RCTD", isample)
    if(opt$marker_genes){
        myRCTD <- readRDS(file.path(RCTD_dir, "results_RCTD_markers.rds"))
    }else{
        myRCTD <- readRDS(file.path(RCTD_dir, "results_RCTD.rds"))
    }
    myRCTD
})

weights <- lapply(RCTD_list, function(myRCTD){
    my_results = myRCTD@results
    my_weights = lapply(my_results, function(x) x$all_weights)
    my_weights_df <- data.frame(do.call(rbind, my_weights))
    my_weights_df
})
weights <- data.frame(do.call(rbind, weights))

coords <- lapply(RCTD_list, function(myRCTD){
    my_coords = myRCTD@spatialRNA@coords
    my_coords
})
coords <- data.frame(do.call(rbind, coords))

spe <- spe[ ,colnames(spe) %in% rownames(coords)]
spe <- spe[ ,match(rownames(coords), colnames(spe))]
colData(spe) <- cbind(colData(spe), weights)

# Add the final clusters
# Add precast results to the spe object
clusters_file <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters", "precast_clusters.csv")
spe[["precast_clusters"]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(clusters_file), by = 'key') |>
    pull(cluster) |>
    as.factor()
# Remove spots with no PRECAST output
spe <- spe[ ,!is.na(spe[["precast_clusters"]])]

# Loop over each sample
for (sample in levels(spe$sample_id)) {
  
  cat("Processing sample:", sample, "\n")
  
  # Subset by sample
  spe_sub <- spe[, spe$sample_id == sample]
  
  # Convert to SingleCellExperiment (required by zellkonverter)
  sce <- as(spe_sub, "SingleCellExperiment")
  
  # Attach spatial coordinates
  coords <- spatialCoords(spe_sub)
  colData(sce)$x <- coords[,2]
  colData(sce)$y <- coords[,1]
  
  # Sample ID, cell types, deconvolution weights, and domain clusters already exist in colData
  # Optional: rename cell type column to "cell_type" if needed
  if (!"cell_type" %in% colnames(colData(sce))) {
    # Assuming RCTD weights are stored with rownames matching cell type names
    # Use max-weight assignment as "cell_type" for convenience
    weights_sub <- colData(sce)[ ,colnames(weights)]
    colData(sce)$cell_type <- colnames(weights_sub)[apply(as.matrix(weights_sub), 1, which.max)]
  }
  
  # Save to .h5ad
  out_file <- file.path(output_dir, paste0(sample, ".h5ad"))
  writeH5AD(sce, out_file)
}