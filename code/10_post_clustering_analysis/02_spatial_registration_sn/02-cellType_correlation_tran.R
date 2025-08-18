library(spatialLIBD)
library(tidyverse)
library(ComplexHeatmap)
library(jaffelab)
library(here)
library(sessioninfo)
library(ggplot2)

opt <- list()
opt$clustering_type <- "01_precast/pseudobulk_capture_area/final_clusters"
#opt$cluster_col <- paste0("BayesSpace_harmony_k", sprintf("%02d", c(3:28)))
opt$cluster_col <- "precast_clusters"
dat_dir <- here("processed-data", "10_post_clustering_analysis", "02_spatial_registration_sn")
plot_dir <- here("plots", "10_post_clustering_analysis", "02_spatial_registration_sn", opt$clustering_type )
spatial_enr_dir <- here("processed-data", "10_post_clustering_analysis","01_pseudobulk_markers",opt$clustering_type)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Load snRNA-seq registration results
sn_registration <- readRDS(file.path(dat_dir, "sn_cellType_registration_tran.rds"))
## Select t-stats from the registration enrichment data
registration_t_stats <- sn_registration$enrichment[, grep("^t_stat", colnames(sn_registration$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))

#### Calculate Correlation Matrix ####
# Get clustering data
modeling_results <- lapply(opt$cluster_col, function(icol){
    get(load(file.path(spatial_enr_dir, paste0("model_results_", icol, ".Rdata"))))
})
names(modeling_results) <- opt$cluster_col

# Compute the correlation between the T-statistics for top 100 genes
cor_top100 <- map(modeling_results, ~ layer_stat_cor(registration_t_stats,
    .x,
    model_type = "enrichment",
    reverse = FALSE,
    top_n = 100
))

if(grepl("precast", opt$clustering_type)){
    if(!grepl("final_clusters", opt$clustering_type)){
        cor_top100 <- lapply(cor_top100, function(icor){
            icor <- icor[ ,paste0("X", sort(as.numeric(gsub("X", "", colnames(icor)))))]
            icor
        })
    }
}

if(grepl("BayesSpace", opt$clustering_type)){
        cor_top100 <- lapply(cor_top100, function(icor){
            icor <- icor[ ,paste0("X", sort(as.numeric(gsub("X", "", colnames(icor)))))]
            icor
        })
}


# Choose colors for visualization
theSeq <- seq(-1, 1, by = 0.01)
my.col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "PRGn"))(length(theSeq))

pdf(file.path(plot_dir, "cor_top100_sn_spatial_registration_final_clusters_tran.pdf"), width = 13, height = 10)
map(cor_top100, ~layer_matrix_plot(t(.x), mypal = my.col, srt = 90, mar = c(10 + max(nchar(rownames(.x)))*0.5, 5 + max(nchar(colnames(.x)))*0.5, 5, 5)), max = 1)
dev.off()


cor_mat <- cor_top100$precast_clusters
rownames(cor_mat) <- gsub("_", " ", rownames(cor_mat))
rownames(cor_mat) <- gsub("\\.", " ", rownames(cor_mat))
colnames(cor_mat) <- gsub("\\.", " ", colnames(cor_mat))
colnames(cor_mat)[colnames(cor_mat) == "Endothelial Ependymal"] <- "Endothelial/\nEpendymal"
order_cols <- c("MSN 1", "MSN 2", "MSN 3", "D1 islands", "Excitatory", "Inhibitory", "WM", "Endothelial/\nEpendymal")
order_rows <- c("MSN D1 A", "MSN D2 A", "MSN D1 D", "MSN D2 B", 
"MSN D1 C", "MSN D2 C", "MSN D2 D", "MSN D1 F", "MSN D1 B", "MSN D1 E", 
"Inhib B", "Inhib A", "Inhib C", "Inhib D", "Inhib E", "Oligo A", "Oligo B",
 "OPC COP", "OPC", "Micro", "Micro resting", "Astro A", "Astro B")

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

pdf(file.path(plot_dir, "cor_top100_sn_spatial_registration_final_clusters_2_tran.pdf"), width = 7, height = 5)
draw(complex_plot_cor, heatmap_legend_side="bottom")
dev.off()

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()