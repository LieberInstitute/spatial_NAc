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
library(ggthemes)

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

spe <- spe[ ,!is.na(spe$precast_final_clusters)]

# Flip some of the samples
spe_Br2743 <- mirrorObject(spe, sample_id = "Br2743", image_id = "lowres", axis = "v")
spe_Br8492 <- mirrorObject(spe, sample_id = "Br8492", image_id = "lowres", axis = "v")
spe_Br8325 <- mirrorObject(spe, sample_id = "Br8325", image_id = "lowres", axis = "v")
spe_Br3942 <- mirrorObject(spe, sample_id = "Br3942", image_id = "lowres", axis = "v")

spe <- spe[ ,!spe$sample_id %in% c("Br2743", "Br8492", "Br8325", "Br3942")]
spe <- cbind(spe, spe_Br2743)
spe <- cbind(spe, spe_Br8492)
spe <- cbind(spe, spe_Br8325)
spe <- cbind(spe, spe_Br3942)
sample_order <- c("Br2743", "Br6432", "Br6423", "Br2720", "Br6471", "Br6522","Br8492", "Br8325", "Br8667", "Br3942")
plot_list <- list()
for (donor in sample_order) {
    plot_list[[donor]] <- spot_plot(
            spe,
            sample_id = donor, var_name = "precast_final_clusters",
            is_discrete = TRUE, spatial = TRUE) + ggtitle(donor)+
            #   Increase size of colored dots in legend
            guides(fill = guide_legend(override.aes = list(size = 5)))
}


pdf(file.path(plot_dir, "clusters_before_merging.pdf"), width = 6.5, height = 6.5)
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

spe$precast_clusters_annotated <- factor(spe$precast_clusters_annotated, 
levels = c("D1 islands", "Endothelial/Ependymal", "Excitatory", "Inhibitory", "MSN 1", "MSN 2", "MSN 3", "WM"))


#safe_colorblind_palette <- c("#652DC1","#0D73B4","#4B5320" ,"#A3EE7D","#EBA31C",  "#3DB1CB", "#D55F0D", "#8F9392")
#safe_colorblind_palette <- palette()
safe_colorblind_palette <- c("#E7298A", "#A6761D","#D95F02" , "#E6AB02",  "#66A61E","#1B9E77", "#7570B3","#666666")
louise_palette <- c("#b2df8a", "#e41a1c", "#377eb8", "#4daf4a", "#ff7f00", "gold", "#a65628",
    "#999999", "black", "grey", "white", "purple")
plot_list <- list()
for (donor in sample_order) {
    plot_list[[donor]] <- spot_plot(
            spe,
            sample_id = donor, var_name = "precast_clusters_annotated",
            is_discrete = TRUE, spatial = TRUE, colors = safe_colorblind_palette) + ggtitle(donor)+
            #   Increase size of colored dots in legend
            guides(fill = guide_legend(override.aes = list(size = 5))) + 
            theme(plot.title = element_text(face = "bold"), legend.position = "right") 
}
pdf(file.path(plot_dir, "clusters_after_merging.pdf"), width = 8, height = 6)
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

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(plot_list[[1]])
p3 <- plot_list[["Br3942"]] + ggtitle("Br3942 (Posterior NAc)") + theme(legend.position = "none")
p2 <- plot_list[["Br6522"]] + ggtitle("Br6522 (Intermediate NAc)")
p1 <- plot_list[["Br6432"]] + ggtitle("Br6432 (Anterior NAc)") + theme(legend.position = "none")

pdf(file.path(plot_dir, "clusters_for_final_figure_Br3942.pdf"), width = 6.8 , height = 6.8)
print(p3)
dev.off()

pdf(file.path(plot_dir, "clusters_for_final_figure_Br6522.pdf"), width = 10, height = 10)
print(p2)
dev.off()

pdf(file.path(plot_dir, "clusters_for_final_figure_Br6432.pdf"), width = 7, height = 7)
print(p1)
dev.off()

# Make plots highlighting individual clusters
clusters_df <- model.matrix(~0+spe$precast_clusters_annotated)
colnames(clusters_df) <- c("D1 islands", "Endothelial/Ependymal", "Excitatory", "Inhibitory", "MSN 1", "MSN 2", "MSN 3", "WM")
colData(spe) <- cbind(colData(spe), data.frame(clusters_df))

colnames(clusters_df) <- gsub(" ", ".", colnames(clusters_df))
colnames(clusters_df) <- gsub("/", ".", colnames(clusters_df))

for(i in colnames(clusters_df)){
    spe[[i]] <- factor(spe[[i]], levels = c(1, 0)) 
    for (donor in sample_order) {
    plot_list[[donor]] <- spot_plot(
            spe,
            sample_id = donor, var_name = i,
            is_discrete = TRUE, spatial = TRUE, colors = c("#006400", "#CCCCCC40")) + ggtitle(donor)+
            #   Increase size of colored dots in legend
            guides(fill = guide_legend(override.aes = list(size = 5))) + 
            theme(plot.title = element_text(face = "bold"), legend.position = "right") 
    }
    pdf(file.path(plot_dir, paste0("clusters_for_final_figure_", i,".pdf")), width = 7, height = 7)
    print(plot_list)
    dev.off()
}