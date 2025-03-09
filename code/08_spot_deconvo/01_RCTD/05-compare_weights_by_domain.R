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

# Read the myRCTD object
spec <- matrix(
    c(
        "marker_genes", "m", 1, "logical", "Use only marker genes from single cell labels?"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

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


# Specify the colors for the domains
safe_colorblind_palette <- c("#66A61E","#1B9E77", "#7570B3","#E7298A","#D95F02" , "#E6AB02","#666666", "#A6761D")
spe$precast_clusters <- factor(spe$precast_clusters, 
levels = c("MSN 1", "MSN 2", "MSN 3", "D1 islands", "Excitatory", "Inhibitory", "WM", "Endothelial/Ependymal"))


plot_df <- data.frame(colData(spe))
plot_df <- plot_df[ ,c("precast_clusters", "DRD1_MSN_A", "DRD2_MSN_A", "DRD2_MSN_B", "DRD1_MSN_C", "DRD1_MSN_B", "DRD1_MSN_D", 
"Excitatory", "Inh_A","Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F", "Astrocyte_A", "Astrocyte_B", "Oligo", "OPC", "Microglia", "Endothelial", "Ependymal")]
plot_df <- reshape2::melt(plot_df, id.vars = "precast_clusters")
colnames(plot_df) <- c("spatial_domain", "cell_type", "RCTD_weights")

plot_df <- plot_df[plot_df$spatial_domain %in% c("MSN 1", "MSN 2", "MSN 3"), ]

plot_df$cell_type <- gsub("_", " ", plot_df$cell_type)
plotDir <- here("plots", "08_spot_deconvo", "01_RCTD")
pdf(file.path(plotDir, "weights_boxplot_marker_genes.pdf"), width = 7, height = 8)
ggplot(plot_df, aes(x = cell_type, y = RCTD_weights, fill = spatial_domain)) + geom_boxplot(outlier.shape = NA) +
coord_flip() + xlab("") +  ylab("RCTD weights") + scale_fill_manual(values = c("#66A61E","#1B9E77", "#7570B3")) +
guides(fill=guide_legend(title="Spatial Domain")) + theme_pubr()
dev.off()

select_cell_types <- c("DRD1 MSN A", "DRD1 MSN B", "DRD1 MSN C", "DRD1 MSN D", "DRD2 MSN A", 
"DRD2 MSN B", "Oligo")
pdf(file.path(plotDir, "weights_boxplot_marker_genes_select_cell_types.pdf"), width = 5, height = 5)
ggplot(plot_df[plot_df$cell_type %in% select_cell_types, ], aes(x = cell_type, y = RCTD_weights, fill = spatial_domain)) + geom_boxplot(outlier.shape = NA) +
coord_flip() + xlab("") +  ylab("RCTD weights") + scale_fill_manual(values = c("#66A61E","#1B9E77", "#7570B3")) +
guides(fill=guide_legend(title="Spatial Domain")) + theme_pubr() + ylim(0, 1)
dev.off()

