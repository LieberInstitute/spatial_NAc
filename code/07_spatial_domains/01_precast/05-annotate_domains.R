library(here)
library(PRECAST)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(spatialNAcUtils)
library(purrr)
library(ggpubr)
library(ggsci)
library(dittoSeq)
library(getopt)
library(pheatmap)

opt <- list()
opt$nnSVG_type <- TRUE
opt$use_random_start <- TRUE
opt$random_start <- 3
opt$K <- 10

if(opt$nnSVG_type){
    if(opt$use_random_start){
        res_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", paste0("random_start_", opt$random_start), "PRECAST_k%s.csv")
        plot_dir <- here("plots", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters")
    }else{
        res_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "PRECAST_k%s.csv")
        plot_dir <- here("plots", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters")
    }
  
}else{
     if(opt$use_random_start){
        res_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_default", paste0("random_start_", opt$random_start), "PRECAST_k%s.csv")
        plot_dir <- here("plots", "07_spatial_domains", "01_precast", "nnSVG_default", "final_clusters")
    }else{
        res_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_default", "PRECAST_k%s.csv")
        plot_dir <- here("plots", "07_spatial_domains", "01_precast", "nnSVG_default", "final_clusters")
    }
}

spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_dimRed_hdf5"
)
spe <- loadHDF5SummarizedExperiment(spe_dir)
# Read in PRECAST results
precast_results <- sprintf(res_path, k = opt$K) |>
        read.csv() |>
        as_tibble() |>
        select(c(key, cluster))

colnames(precast_results) <- c("key", "precast_final_clusters")

temp <- colnames(spe)
colData(spe) <- colData(spe) |>
    as_tibble() |>
    left_join(precast_results, by = "key") |>
    DataFrame()
colnames(spe) <- temp

plot_list <- list()
for (donor in levels(spe$sample_id)) {
    plot_list[[donor]] <- spot_plot(
            spe,
            sample_id = donor, var_name = "precast_final_clusters",
            is_discrete = TRUE, spatial = TRUE) + ggtitle(donor)+
            #   Increase size of colored dots in legend
            guides(fill = guide_legend(override.aes = list(size = 5)))
}

pdf(file.path(plot_dir, "clusters_before_merging.pdf"), width = 5, height = 5)
print(plot_list)
dev.off()

# Merge and annotate clusters
spe$precast_clusters_annotated <- spe$precast_final_clusters
spe$precast_clusters_annotated[spe$precast_clusters_annotated == 1] <- "MSN 1"
spe$precast_clusters_annotated[spe$precast_clusters_annotated == 2] <- "WM"
spe$precast_clusters_annotated[spe$precast_clusters_annotated == 3] <- "D1 islands"
spe$precast_clusters_annotated[spe$precast_clusters_annotated == 4] <- "MSN 2"
spe$precast_clusters_annotated[spe$precast_clusters_annotated == 5] <- "MSN 3"
spe$precast_clusters_annotated[spe$precast_clusters_annotated == 6] <- "Endothelial/Ependymal"
spe$precast_clusters_annotated[spe$precast_clusters_annotated == 7] <- "Excitatory"
spe$precast_clusters_annotated[spe$precast_clusters_annotated == 8] <- "Inhibitory"
spe$precast_clusters_annotated[spe$precast_clusters_annotated == 9] <- "WM"
spe$precast_clusters_annotated[spe$precast_clusters_annotated == 10] <- "Endothelial/Ependymal"

spe$precast_clusters_annotated <- factor(spe$precast_clusters_annotated, levels = c("MSN 1", "MSN 2", "MSN 3", "D1 islands", "Excitatory", "Inhibitory", "Endothelial/Ependymal", "WM"))

plot_list <- list()
for (donor in levels(spe$sample_id)) {
    plot_list[[donor]] <- spot_plot(
            spe,
            sample_id = donor, var_name = "precast_clusters_annotated",
            is_discrete = TRUE, spatial = TRUE) + ggtitle(donor)+
            #   Increase size of colored dots in legend
            guides(fill = guide_legend(override.aes = list(size = 5)))
}
pdf(file.path(plot_dir, "clusters_after_merging.pdf"), width = 10, height = 10)
print(plot_list)
dev.off()

clusters_csv <- colData(spe)[ ,c("key", "precast_clusters_annotated")]
if(opt$nnSVG_type){
    saveDir <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters")
}else{
    saveDir <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_default", "final_clusters")
}
dir.create(saveDir, recursive = TRUE, showWarnings = FALSE)
colnames(clusters_csv) <- c("key", "cluster")
write.csv(clusters_csv, file.path(saveDir, "precast_clusters.csv"), row.names = FALSE, quote = FALSE)