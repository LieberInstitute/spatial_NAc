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
codeDir <- here("code")
source(file.path(codeDir, "plot_utils.R"))

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

spe <- spe[ ,colnames(spe) %in% rownames(coords)]
spe <- spe[ ,match(rownames(coords), colnames(spe))]
colData(spe) <- cbind(colData(spe), weights_df)

celltypes = levels(myRCTD@reference@cell_types)
plot_list <- lapply(celltypes, function(i) {
  print(i)  
  vis_gene(spe, sampleid = sample_id, geneid = i, is_stitched = TRUE, cont_colors = viridisLite::rocket(10, direction = -1))
})

if(opt$marker_genes){
    plot_dir <- here("plots", "08_spot_deconvo", "01_RCTD", sample_id, "marker_genes")
}else{
    plot_dir <- here("plots", "08_spot_deconvo", "01_RCTD", sample_id, "all_genes")
}

plot_list[[1]] <- plot_list[[1]] + ggtitle("Oligodendrocytes") + theme(plot.title = element_text(face = "bold"))
plot_list[[2]] <- plot_list[[2]] + ggtitle("DRD1 MSN A") + theme(plot.title = element_text(face = "bold"))
plot_list[[3]] <- plot_list[[3]] + ggtitle("Microglia") + theme(plot.title = element_text(face = "bold"))
plot_list[[4]] <- plot_list[[4]] + ggtitle("OPC") + theme(plot.title = element_text(face = "bold"))
plot_list[[5]] <- plot_list[[5]] + ggtitle("DRD2 MSN A") + theme(plot.title = element_text(face = "bold"))
plot_list[[6]] <- plot_list[[6]] + ggtitle("Ependymal") + theme(plot.title = element_text(face = "bold"))
plot_list[[7]] <- plot_list[[7]] + ggtitle("Astrocyte A") + theme(plot.title = element_text(face = "bold"))
plot_list[[8]] <- plot_list[[8]] + ggtitle("DRD1 MSN B") + theme(plot.title = element_text(face = "bold"))
plot_list[[9]] <- plot_list[[9]] + ggtitle("Endothelial") + theme(plot.title = element_text(face = "bold"))
plot_list[[10]] <- plot_list[[10]] + ggtitle("Inhibitory A") + theme(plot.title = element_text(face = "bold"))
plot_list[[11]] <- plot_list[[11]] + ggtitle("DRD2 MSN B") + theme(plot.title = element_text(face = "bold"))
plot_list[[12]] <- plot_list[[12]] + ggtitle("Astrocyte B") + theme(plot.title = element_text(face = "bold"))
plot_list[[13]] <- plot_list[[13]] + ggtitle("DRD1 MSN C") + theme(plot.title = element_text(face = "bold"))
plot_list[[14]] <- plot_list[[14]] + ggtitle("DRD1 MSN D") + theme(plot.title = element_text(face = "bold"))
plot_list[[15]] <- plot_list[[15]] + ggtitle("Inhibitory B") + theme(plot.title = element_text(face = "bold"))
plot_list[[16]] <- plot_list[[16]] + ggtitle("Inhibitory C") + theme(plot.title = element_text(face = "bold"))
plot_list[[17]] <- plot_list[[17]] + ggtitle("Inhibitory D") + theme(plot.title = element_text(face = "bold"))
plot_list[[18]] <- plot_list[[18]] + ggtitle("Inhibitory E") + theme(plot.title = element_text(face = "bold"))
plot_list[[19]] <- plot_list[[19]] + ggtitle("Excitatory") + theme(plot.title = element_text(face = "bold"))
plot_list[[20]] <- plot_list[[20]] + ggtitle("Inhibitory F") + theme(plot.title = element_text(face = "bold"))


dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(plot_dir,"multi_allWeight.pdf"), width = 7, height = 7)
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

colnames(df) <- paste0(colnames(df), "_2")
colData(spe) <- cbind(colData(spe), df)

plot_list_2 <- lapply(celltypes, function(i) {
    print(i)  
    vis_gene(spe, sampleid = sample_id, geneid = paste0(i, "_2"), is_stitched = TRUE, cont_colors = viridisLite::rocket(10, direction = -1))
})

plot_list_2[[1]] <- plot_list_2[[1]] + ggtitle("Oligodendrocytes") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[2]] <- plot_list_2[[2]] + ggtitle("DRD1 MSN A") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[3]] <- plot_list_2[[3]] + ggtitle("Microglia") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[4]] <- plot_list_2[[4]] + ggtitle("OPC") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[5]] <- plot_list_2[[5]] + ggtitle("DRD2 MSN A") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[6]] <- plot_list_2[[6]] + ggtitle("Ependymal") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[7]] <- plot_list_2[[7]] + ggtitle("Astrocyte A") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[8]] <- plot_list_2[[8]] + ggtitle("DRD1 MSN B") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[9]] <- plot_list_2[[9]] + ggtitle("Endothelial") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[10]] <- plot_list_2[[10]] + ggtitle("Inhibitory A") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[11]] <- plot_list_2[[11]] + ggtitle("DRD2 MSN B") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[12]] <- plot_list_2[[12]] + ggtitle("Astrocyte B") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[13]] <- plot_list_2[[13]] + ggtitle("DRD1 MSN C") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[14]] <- plot_list_2[[14]] + ggtitle("DRD1 MSN D") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[15]] <- plot_list_2[[15]] + ggtitle("Inhibitory B") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[16]] <- plot_list_2[[16]] + ggtitle("Inhibitory C") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[17]] <- plot_list_2[[17]] + ggtitle("Inhibitory D") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[18]] <- plot_list_2[[18]] + ggtitle("Inhibitory E") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[19]] <- plot_list_2[[19]] + ggtitle("Excitatory") + theme(plot.title = element_text(face = "bold"))
plot_list_2[[20]] <- plot_list_2[[20]] + ggtitle("Inhibitory F") + theme(plot.title = element_text(face = "bold"))


pdf(file.path(plot_dir,"multi_subWeight_markers.pdf"), width = 6, height = 6)
for(i in c(1:length(plot_list_2))){
    print(plot_list_2[[i]])
}
dev.off()