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
library(viridis)
library(ggpubr)
library(tidyverse)

opt <- list()
opt$marker_genes <- TRUE
spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5")
spe <- loadHDF5SummarizedExperiment(spe_dir)
sample_ids <- levels(spe$sample_id)

# Read the RCTD results
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

# Visualize the composition of the D1 islands
plot_df <- data.frame(colData(spe))
plot_df <- plot_df[ ,c("precast_clusters", "DRD1_MSN_A", "DRD2_MSN_A", "DRD2_MSN_B", "DRD1_MSN_C", "DRD1_MSN_B", "DRD1_MSN_D", 
"Excitatory", "Inh_A","Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F", "Astrocyte_A", "Astrocyte_B", "Oligo", "OPC", "Microglia", "Endothelial", "Ependymal")]
plot_df <- reshape2::melt(plot_df, id.vars = "precast_clusters")
colnames(plot_df) <- c("spatial_domain", "cell_type", "RCTD_weights")

plot_df <- plot_df[plot_df$spatial_domain %in% c("D1 islands"), ]
plot_df$cell_type <- gsub("_", " ", plot_df$cell_type)

mean_by_group <- plot_df %>%
  group_by(cell_type) %>%
  summarize(mean_value = mean(RCTD_weights))

select_cell_types <- mean_by_group$cell_type[mean_by_group$mean_value > 0.06]

plotDir <- here("plots", "08_spot_deconvo", "01_RCTD")
pdf(file.path(plotDir, "RCTD_weights_boxplot_D1_islands.pdf"), width = 4, height = 4)
ggplot(plot_df[plot_df$cell_type %in% select_cell_types, ], aes(x = cell_type, y = RCTD_weights, fill = cell_type)) + geom_boxplot(outlier.shape = NA) +
coord_flip() + xlab("") +  ylab("RCTD weights") + theme_pubr() + theme(legend.position = "none") 
dev.off()

plot_df <- data.frame(colData(spe))
plot_df <- plot_df[ ,c("precast_clusters", "sample_id", "DRD1_MSN_A", "DRD2_MSN_A", "DRD2_MSN_B", "DRD1_MSN_C", "DRD1_MSN_B", "DRD1_MSN_D", 
"Excitatory", "Inh_A","Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F", "Astrocyte_A", "Astrocyte_B", "Oligo", "OPC", "Microglia", "Endothelial", "Ependymal")]
plot_df <- plot_df[plot_df$precast_clusters == "D1 islands", ]
plot_df <- plot_df[ ,!colnames(plot_df) == "precast_clusters"]
plot_df <- reshape2::melt(plot_df, id.vars = "sample_id")
plot_df <- plot_df[plot_df$sample_id %in% c("Br2720", "Br6471", "Br6522", "Br8492"), ]
colnames(plot_df) <- c("sample_id", "cell_type", "RCTD_weights")
plot_df$cell_type <- gsub("_", " ", plot_df$cell_type)

plot_df <- plot_df[plot_df$cell_type %in% c("DRD1 MSN B", "DRD1 MSN C", "DRD1 MSN D"), ]
pdf(file.path(plotDir, "D1_islands_weights_by_sample.pdf"), width = 6, height = 5)
ggplot(plot_df, aes(x = sample_id, y = RCTD_weights, fill = cell_type)) + geom_boxplot(outlier.shape = NA) +
coord_flip() + xlab("") +  ylab("RCTD weights") + theme_pubr()
dev.off()