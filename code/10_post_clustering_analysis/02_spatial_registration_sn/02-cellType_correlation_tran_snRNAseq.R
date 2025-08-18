library(spatialLIBD)
library(tidyverse)
library(ComplexHeatmap)
library(jaffelab)
library(here)
library(sessioninfo)
library(ggplot2)

dat_dir <- here("processed-data", "10_post_clustering_analysis", "02_spatial_registration_sn")
plot_dir <- here("plots", "10_post_clustering_analysis", "02_spatial_registration_sn")
# Load snRNA-seq registration results
subset_neurons <- TRUE
if(subset_neurons){
    sn_registration <- readRDS(file.path(dat_dir, "sn_cellType_registration_tran_neurons.rds"))
}else{
    sn_registration <- readRDS(file.path(dat_dir, "sn_cellType_registration_tran.rds"))
}

## Select t-stats from the registration enrichment data
registration_t_stats <- sn_registration$enrichment[, grep("^t_stat", colnames(sn_registration$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))

if(subset_neurons){
    sn_file_name <- c("sn_cellType_registration_neurons.rds")
}else{
    sn_file_name <- c("sn_cellType_registration.rds")
}


#### Calculate Correlation Matrix ####
# Get clustering data
modeling_results <- lapply(sn_file_name, function(ifile){
    readRDS(file.path(dat_dir, ifile))
})
names(modeling_results) <- sn_file_name

# Compute the correlation between the T-statistics for top 100 genes
cor_top100 <- map(modeling_results, ~ layer_stat_cor(registration_t_stats,
    .x,
    model_type = "enrichment",
    reverse = FALSE,
    top_n = 100
))


# Choose colors for visualization
theSeq <- seq(-1, 1, by = 0.01)
my.col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "PRGn"))(length(theSeq))

if(subset_neurons){
    pdf(file.path(plot_dir, "cor_top100_sn_final_clusters_tran_neurons.pdf"), width = 13, height = 10)
    map(cor_top100, ~layer_matrix_plot(t(.x), mypal = my.col, srt = 90, mar = c(10 + max(nchar(rownames(.x)))*0.5, 5 + max(nchar(colnames(.x)))*0.5, 5, 5)), max = 1)
    dev.off()
}else{
    pdf(file.path(plot_dir, "cor_top100_sn_final_clusters_tran.pdf"), width = 13, height = 10)
    map(cor_top100, ~layer_matrix_plot(t(.x), mypal = my.col, srt = 90, mar = c(10 + max(nchar(rownames(.x)))*0.5, 5 + max(nchar(colnames(.x)))*0.5, 5, 5)), max = 1)
    dev.off()
}

if(subset_neurons){
    cor_mat <- cor_top100$sn_cellType_registration_neurons.rds
    rownames(cor_mat) <- gsub("_", " ", rownames(cor_mat))
    rownames(cor_mat) <- gsub("\\.", " ", rownames(cor_mat))
    colnames(cor_mat) <- gsub("_", " ", colnames(cor_mat))
    order_cols <- c("DRD1 MSN A", "DRD2 MSN A", "DRD1 MSN C", "DRD2 MSN B", "DRD1 MSN B", "DRD1 MSN D", "Excitatory", "Inh A", 
    "Inh B", "Inh C", "Inh D", "Inh E", "Inh F")
    order_rows <- c("MSN D1 A", "MSN D2 A", "MSN D1 D", "MSN D2 B", 
    "MSN D1 C","MSN D2 D", "MSN D2 C", "MSN D1 E", "MSN D1 B", "MSN D1 F", "Inhib C", "Inhib D", 
    "Inhib B", "Inhib A", "Inhib E")

    cor_mat <- cor_mat[rownames(cor_mat) %in% order_rows, ]
    cor_mat <- cor_mat[match(order_rows, rownames(cor_mat)), ]
    cor_mat <- cor_mat[ ,match(order_cols, colnames(cor_mat))]

    col_fun <- circlize::colorRamp2(c(-1.0,0,1.0),hcl_palette = "Purple-Green")
    complex_plot_cor <- ComplexHeatmap::Heatmap(matrix = t(cor_mat),
                              name = "Correlation",
                              column_title = NULL,
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              row_title = NULL,
                              rect_gp = gpar(col = "white", lwd = 2),
                              col = col_fun, 
                              border = TRUE,
                              heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(6, "cm"),title_position = "topcenter", title_gp = gpar(fontsize = 14), border = "black"), 
                              row_names_side = "left")

    pdf(file.path(plot_dir, "cor_top100_sn_final_clusters2_tran_neurons.pdf"), width = 7, height = 6)
    draw(complex_plot_cor, heatmap_legend_side="bottom")
    dev.off()

}else{
    cor_mat <- cor_top100$sn_cellType_registration.rds
    rownames(cor_mat) <- gsub("_", " ", rownames(cor_mat))
    rownames(cor_mat) <- gsub("\\.", " ", rownames(cor_mat))
    colnames(cor_mat) <- gsub("_", " ", colnames(cor_mat))
    order_cols <- c("DRD1 MSN A", "DRD2 MSN A", "DRD1 MSN C", "DRD2 MSN B", "DRD1 MSN B", "DRD1 MSN D", "Excitatory", "Inh A", 
    "Inh B", "Inh C", "Inh D", "Inh E", "Inh F", "Astrocyte A", "Astrocyte B", "Microglia", "OPC", "Oligo", "Ependymal", "Endothelial")
    order_rows <- c("MSN D1 A", "MSN D2 A", "MSN D1 D", "MSN D2 B", 
    "MSN D1 C","MSN D2 D", "MSN D2 C", "MSN D1 E", "MSN D1 B", "MSN D1 F", "Inhib C", "Inhib D", 
    "Inhib B", "Inhib A", "Inhib E","Astro B", "Astro A",
     "Micro", "Micro resting", "OPC","OPC COP",
    "Oligo A", "Oligo B")


    cor_mat <- cor_mat[rownames(cor_mat) %in% order_rows, ]
    cor_mat <- cor_mat[match(order_rows, rownames(cor_mat)), ]
    cor_mat <- cor_mat[ ,match(order_cols, colnames(cor_mat))]

    col_fun <- circlize::colorRamp2(c(-1.0,0,1.0),hcl_palette = "Purple-Green")
    complex_plot_cor <- ComplexHeatmap::Heatmap(matrix = t(cor_mat),
                              name = "Correlation",
                              column_title = NULL,
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              row_title = NULL,
                              rect_gp = gpar(col = "white", lwd = 2),
                              col = col_fun, 
                              border = TRUE,
                              heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(6, "cm"),title_position = "topcenter", title_gp = gpar(fontsize = 14), border = "black"), 
                              row_names_side = "left")

    pdf(file.path(plot_dir, "cor_top100_sn_final_clusters2_tran.pdf"), width = 7, height = 6)
    draw(complex_plot_cor, heatmap_legend_side="bottom")
    dev.off()
}


print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()