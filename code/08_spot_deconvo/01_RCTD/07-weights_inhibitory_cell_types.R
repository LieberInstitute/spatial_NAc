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

# Visualize the distribution of inhibitory cell types
plot_df <- data.frame(colData(spe))
plot_df <- plot_df[ ,c("precast_clusters", "Inh_A","Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F")]
plot_df <- reshape2::melt(plot_df, id.vars = "precast_clusters")
colnames(plot_df) <- c("spatial_domain", "cell_type", "RCTD_weights")

plotDir <- here("plots", "08_spot_deconvo", "01_RCTD")

filtered_df <- plot_df %>% filter(RCTD_weights > 0.1)
safe_colorblind_palette <- c(
  "#E7298A", "#A6761D", "#D95F02", "#E6AB02",
  "#66A61E", "#1B9E77", "#7570B3", "#666666"
)


pdf(file.path(plotDir, "RCTD_weights_Inhib_cellTypes.pdf"), width = 12, height = 8)
ggplot(filtered_df, aes(x = spatial_domain, y = RCTD_weights, fill = spatial_domain)) +
  geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.95) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5, color = "black") +
  facet_wrap(~ cell_type, scales = "free_y") +
  scale_fill_manual(values = safe_colorblind_palette) +
  labs(
    title = "Distribution of RCTD Weights > 0.1 by Inhibitory Cell Type and Spatial Domain",
    x = "Spatial Domain",
    y = "RCTD Weight"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    panel.grid = element_blank(),  # remove gridlines
    panel.background = element_rect(fill = "white", color = "black"),  # individual plot boxes
    strip.background = element_rect(fill = "white", color = "black"),
    legend.position = "none"
  )
dev.off()

inh_cellTypes <- c("Inh_A", "Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F")
sample_order <- c("Br2743", "Br6432", "Br6423", "Br2720", "Br6471", "Br6522","Br8492", "Br8325", "Br8667", "Br3942")
spe$precast_clusters <- factor(as.character(spe$precast_clusters), levels = c("D1 islands", "Endothelial/Ependymal", "Excitatory", "Inhibitory", "MSN 1", "MSN 2", "MSN 3", "WM")) 
safe_colorblind_palette <- c("#E7298A", "#A6761D","#D95F02" , "#E6AB02",  "#66A61E","#1B9E77", "#7570B3","#666666")
names(safe_colorblind_palette) <- levels(spe$precast_clusters)
for (cellType in inh_cellTypes) {
  print(cellType)
  plot_list <- list()
  for (donor in sample_order) {
    spe_sub <- spe[ ,spe$sample_id == donor]
    spe_sub <- spe_sub[ ,!spe_sub$exclude_overlapping]
    plot_list[[donor]] <- make_escheR(spe_sub) |> 
                          add_fill(var = cellType) |> 
                          add_ground(var = "precast_clusters", stroke = 0.5) +
                          scale_color_manual(values = safe_colorblind_palette) +
                          scale_fill_gradient(low = "white", high = "black")
  }
  pdf(file.path(plotDir, paste0(cellType, "_spot_plots_escheR.pdf")), width = 12, height = 12)
  for (p in plot_list) print(p)
  dev.off()
}

plot_list <- list()
donor <- "Br6522"
for (cellType in inh_cellTypes) {
    spe_sub <- spe[ ,spe$sample_id == donor]
    spe_sub <- spe_sub[ ,!spe_sub$exclude_overlapping]
    plot_list[[cellType]] <- make_escheR(spe_sub) |> 
                          add_fill(var = cellType) |> 
                          add_ground(var = "precast_clusters", stroke = 0.5) +
                          scale_color_manual(values = safe_colorblind_palette) +
                          scale_fill_gradient(low = "white", high = "black")
}
pdf(file.path(plotDir, "Br6522_Inhib_spot_plots_escheR.pdf"), width = 11, height = 11)
for (p in plot_list) print(p)
dev.off()
