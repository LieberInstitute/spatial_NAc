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
spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5")
spe <- loadHDF5SummarizedExperiment(spe_dir)
sample_id <- levels(spe$sample_id)[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]

# Read in the myRCTD object
RCTD_dir <- here::here("processed-data", "08_spot_deconvo", "01_RCTD", sample_id)
if(opt$marker_genes){
    myRCTD <- readRDS(file.path(RCTD_dir, "results_RCTD_markers.rds"))
}else{
    myRCTD <- readRDS(file.path(RCTD_dir, "results_RCTD.rds"))
}

## multi mode ## 
print('Examining multi mode all weights results')

results = myRCTD@results
weights = lapply(results, function(x) x$all_weights)
weights_df <- data.frame(do.call(rbind, weights))
norm_weights <- normalize_weights(weights_df)
coords = myRCTD@spatialRNA@coords

celltypes = levels(myRCTD@reference@cell_types)
plot_list <- lapply(celltypes, function(i) {
  print(i)  
  ggplot(coords, aes(x = coords$x, y=coords$y, color = weights_df[,i])) + labs(title = i, x="", y="") + 
    geom_point(size = 0.5)+scale_color_gradientn(colours = viridis(10, option = "magma"), limits = c(0,1)) +
    scale_y_reverse()+ theme_pubr() + theme(legend.key.width = unit(1, "cm"))+labs(color = "") 
})

if(opt$marker_genes){
    plot_dir <- here("plots", "08_spot_deconvo", "01_RCTD", sample_id, "marker_genes")
}else{
    plot_dir <- here("plots", "08_spot_deconvo", "01_RCTD", sample_id, "all_genes")
}

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(plot_dir,"multi_allWeight.pdf"), width = 6, height = 6)
for(i in c(1:length(plot_list))){
    print(plot_list[[i]])
}
dev.off()

weights_df = as.data.frame(weights_df)
rownames(weights_df) = rownames(myRCTD@spatialRNA@coords)
if(opt$marker_genes){
    write.csv(weights_df, file.path(RCTD_dir, "clusters_allWeights_markers.csv"))
}else{
    write.csv(weights_df, file.path(RCTD_dir, "clusters_allWeights.csv"))
}


print('Examining multi mode sub weights results')
results = myRCTD@results
weights = lapply(results, function(x) x$sub_weights)
df = weights[[1]]
missing_columns <- setdiff(celltypes, names(df))
if (length(missing_columns) > 0) {for (col_name in missing_columns) {df[[col_name]] <- 0 }}
df <- as.data.frame(t(df))
if (length(missing_columns) == length(celltypes)) {df = select(df, -c(setdiff(colnames(df), celltypes)))}

for (i in 2:length(weights)) {
  # Extract the list
  current_list <- weights[[i]]
  missing_columns <- setdiff(celltypes, names(current_list))
  if (length(missing_columns) > 0) {for (col_name in missing_columns) {current_list[[col_name]] <- 0 }}
  current_list = as.data.frame(t(current_list))
  if (length(missing_columns) == length(celltypes)) {current_list = select(current_list, -c(setdiff(colnames(current_list), celltypes)))}
  current_list = current_list[,colnames(df)]
  # Combine the current list with the combined dataframe by column names
  df <- rbind(df, current_list)
}

plot_list <- lapply(celltypes, function(i) {
  ggplot(coords, aes(x = coords$x, y=coords$y , color = df[,i])) + labs(title = i, x="", y="") + 
    geom_point(size = 0.5)+scale_color_gradientn(colours = viridis(10, option = "magma"), limits = c(0,1)) +
    scale_y_reverse() + theme_pubr() + theme(legend.key.width = unit(1, "cm"))+labs(color = "")
})

pdf(file.path(plot_dir,"multi_subWeight_markers.pdf"), width = 6, height = 6)
for(i in c(1:length(plot_list))){
    print(plot_list[[i]])
}
dev.off()