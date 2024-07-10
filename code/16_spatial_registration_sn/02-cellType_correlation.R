library(spatialLIBD)
library(tidyverse)
library(ComplexHeatmap)
library(jaffelab)
library(here)
library(sessioninfo)
library(ggplot2)

dat_dir <- here("processed-data", "16_spatial_registration_sn")
plot_dir <- here("plots", "16_spatial_registration_sn")
spatial_enr_dir <- here("processed-data", "14_pseudobulk_spatial", "01_precast", "pseudobulk_capture_area")

# Load snRNA-seq registration results
sn_registration <- readRDS(file.path(dat_dir, "sn_cellType_registration.rds"))
## Select t-stats from the registration enrichment data
registration_t_stats <- sn_registration$enrichment[, grep("^t_stat", colnames(sn_registration$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))

#### Calculate Correlation Matrix ####
# Get clustering data
modeling_results <- lapply(c(3:18), function(K){
    get(load(file.path(spatial_enr_dir, paste0("model_results_precast_k", K, ".Rdata"))))
})
names(modeling_results) <- paste0("K_", c(3:8))

# Compute the correlation between the T-statistics for top 100 genes
cor_top100 <- map(modeling_results, ~ layer_stat_cor(registration_t_stats,
    .x,
    model_type = "enrichment",
    reverse = FALSE,
    top_n = 100
))

cor_top100 <- lapply(cor_top100, function(icor){
    icor <- icor[ ,paste0("X", sort(as.numeric(gsub("X", "", colnames(icor)))))]
    icor
})

# Choose colors for visualization
theSeq <- seq(-1, 1, by = 0.01)
my.col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "PRGn"))(length(theSeq))

pdf(file.path(plot_dir, "cor_top100_sn_spatial_registration.pdf"), width = 12, height = 10)
map(cor_top100, ~layer_matrix_plot(t(.x), mypal = my.col, srt = 90, mar = c(10 + max(nchar(rownames(.x)))*0.5, 5, 5, 5)), max = 1)
dev.off()


## Annotate
layer_anno_easy <- map2(cor_top100, names(cor_top100), function(cor, name) {
    anno <- annotate_registered_clusters(
        cor_stats_layer = cor,
        confidence_threshold = 0.25,
        cutoff_merge_ratio = 0.25
    )
    colnames(anno) <- gsub("layer", name, colnames(anno))
    return(anno)
})

layer_anno_strict <- map2(cor_top100, names(cor_top100), function(cor, name) {
    anno <- annotate_registered_clusters(
        cor_stats_layer = cor,
        confidence_threshold = 0.25,
        cutoff_merge_ratio = 0.1
    )
    colnames(anno) <- gsub("layer", name, colnames(anno))
    return(anno)
})

# Differences in layer annotations between easy and strict
layer_anno_easy$K_10 |>
    left_join(layer_anno_strict$K_10 |> select(cluster, strict_label = K_10_label)) |>
    filter(K_10_label != strict_label)


