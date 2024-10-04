library(here)
library(PRECAST)
library(HDF5Array)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(sessioninfo)
library(tidyverse)
library(Matrix)
library(SpatialExperiment)
library(spatialNAcUtils)
library(scran)
library(scater)
# Fixed number of clusters aimed at identifying WM and GM
k <- 2

spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5"
)
hvg_path <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "top.hvgs.Rdata")
out_path <- here("processed-data", "10_precast/00_pre_clustering", paste0("PRECAST_k", k, ".csv"))
plot_dir <- here('plots', '10_precast', '00_pre_clustering')

set.seed(1)
dir.create(dirname(out_path), showWarnings = FALSE)
dir.create(dirname(plot_dir), showWarnings = FALSE)
spe <- loadHDF5SummarizedExperiment(spe_dir)

#   PRECAST expects array coordinates in 'row' and 'col' columns
rownames(spe) <- rowData(spe)$gene_name
rownames(spe) <- make.names(rownames(spe), unique = TRUE)
spe$row <- spe$array_row_transformed
spe$col <- spe$array_col_transformed

#   Create a list of Seurat objects: one per donor
seu_list <- lapply(
    levels(spe$donor),
    function(donor) {
        cat(donor, "\n")
        small_spe <- spe[, spe$donor == donor]
        cat(dim(small_spe), "\n")
        CreateSeuratObject(
            #   Bring into memory to greatly improve speed
            counts = as(assays(small_spe)$counts, "dgCMatrix"),
            meta.data = as.data.frame(colData(small_spe)),
            project = "spatialNAc"
        )
    }
)

load(hvg_path)
hvgs_id <- top.hvgs.p2
hvgs_name <- c()
for(i in c(1:length(top.hvgs.p2))){
    hvgs_name[i] <- rowData(spe)$gene_name[rowData(spe)$gene_id == top.hvgs.p2[i]]
}

pre_obj <- CreatePRECASTObject(
    seuList = seu_list,
    selectGenesMethod = NULL,
    customGenelist = hvgs_name)
    #   Using defaults for gene-filtering-related parameters. Though each donor
    #   consists of more spots than 1 typical Visium capture area (and would
    #   thus be expected to throw off the appropriateness of the defaults for
    #   'premin.spots', etc), we're using top 20% of the HVGs
#   Setting platform to "Visium" just means to use array indices, which should
#   work fine despite the abnormal/ "artificial" capture area we've created by
#   stitching
pre_obj <- AddAdjList(pre_obj, platform = "Visium")

#   Following https://feiyoung.github.io/PRECAST/articles/PRECAST.BreastCancer.html,
#   which involves overriding some default values, though the implications are not
#   documented
pre_obj <- AddParSetting(
    pre_obj,
    Sigma_equal = FALSE, verbose = TRUE, maxIter = 30
)

#   Fit model
pre_obj <- PRECAST(pre_obj, K = k)
## backup the fitting results in resList
resList <- PRECASTObj@resList
pre_obj <- SelectModel(pre_obj)
pre_obj <- IntegrateSpaData(pre_obj, species = "Human")

cols_cluster <- chooseColors(palettes_name = "Classic 20", n_colors = 7, plot_colors = TRUE)[c(1, 3)]

pre_obj <- RunTSNE(pre_obj, reduction = "PRECAST", tSNE.method = "FIt-SNE")

p1 <- dimPlot(pre_obj, item = "cluster", point_size = 0.5, font_family = "serif", cols = cols_cluster,
    border_col = "gray10", nrow.legend = 14, legend_pos = "right")

cols_batch <- chooseColors(palettes_name = "Classic 20", n_colors = 10, plot_colors = TRUE)
p2 <- dimPlot(pre_obj, item = "batch", point_size = 0.5, font_family = "serif", cols = cols_batch,
    border_col = "gray10", nrow.legend = 14, legend_pos = "right")

pdf(file.path(plot_dir, "/precast_tSNE.pdf"), width = 12, height = 8)
plot_grid(p1, p2, ncol = 2)
dev.off()

dat_deg <- FindAllMarkers(pre_obj)
n <- 20
dat_deg %>%
    group_by(cluster) %>%
    top_n(n = n, wt = avg_log2FC) -> top40

pdf(file.path(plot_dir, "/precast_cluster_markers.pdf"), width = 15, height = 7)
DotPlot(pre_obj, features = c(unique(top40$gene), "MBP", "MOBP"), col.min = 0, col.max = 1) + theme(axis.text.x = element_text(angle = 45,
    hjust = 1, size = 8))
dev.off()

#   Extract PRECAST results, clean up column names, and export to CSV

precast_results <- pre_obj@meta.data |>
    rownames_to_column("key") |>
    as_tibble() |>
    select(-orig.ident) |>
    rename_with(~ sub("_PRE_CAST", "", .x))
    
write_csv(precast_results, out_path)

precast_cluster <- colData(spe) |>
        as_tibble() |>
        left_join(precast_results, by = "key") |>
        pull(cluster)

spe$precast_cluster <- precast_cluster
# Extract PRECAST results 
pdf(file.path(plot_dir, "WM_clusters.pdf"), width = 8, height = 8)
spot_plot(spe, "Br2720", var_name = "precast_cluster", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br2743", var_name = "precast_cluster", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br3942", var_name = "precast_cluster", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6423", var_name = "precast_cluster", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6432", var_name = "precast_cluster", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6471", var_name = "precast_cluster", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6522", var_name = "precast_cluster", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br8325", var_name = "precast_cluster", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br8492", var_name = "precast_cluster", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br8667", var_name = "precast_cluster", is_discrete = TRUE, spatial = TRUE)
dev.off()

df <- data.frame(colData(spe))
Outlier_df <- df[df$low_umi, ]
nonOutlier_df <- df[!df$low_umi, ]

Outlier_summ <- Outlier_df %>% group_by(sample_id, precast_cluster, exclude_overlapping) %>% tally()
nonOutlier_summ <- nonOutlier_df %>% group_by(sample_id, precast_cluster, exclude_overlapping) %>% tally()

pdf(file.path(plot_dir, "Outlier_Low_UMI_spot_distribution.pdf"), width = 10, height = 6)
ggplot(Outlier_summ, aes(x = sample_id, y = n, fill = exclude_overlapping)) + geom_bar(stat="identity", position = position_dodge(preserve = "single")) +
facet_wrap(~precast_cluster, nrow = 3) + theme_classic() + ggtitle("Outlier (Low UMI) spot distribution")
dev.off()

pdf(file.path(plot_dir, "nonOutlier_Low_UMI_spot_distribution.pdf"), width = 10, height = 6)
ggplot(nonOutlier_summ, aes(x = sample_id, y = n, fill = exclude_overlapping)) + geom_bar(stat="identity", position = position_dodge(preserve = "single")) +
facet_wrap(~precast_cluster, nrow = 3) + theme_classic() + ggtitle("Non-outlier (Low UMI) spot distribution")
dev.off()

Outlier_df <- df[df$low_gene, ]
nonOutlier_df <- df[!df$low_gene, ]

Outlier_summ <- Outlier_df %>% group_by(sample_id, precast_cluster, exclude_overlapping) %>% tally()
nonOutlier_summ <- nonOutlier_df %>% group_by(sample_id, precast_cluster, exclude_overlapping) %>% tally()

pdf(file.path(plot_dir, "Outlier_Low_gene_spot_distribution.pdf"), width = 10, height = 6)
ggplot(Outlier_summ, aes(x = sample_id, y = n, fill = exclude_overlapping)) + geom_bar(stat="identity", position = position_dodge(preserve = "single")) +
facet_wrap(~precast_cluster, nrow = 3) + theme_classic() + ggtitle("Outlier (Low gene) spot distribution")
dev.off()

pdf(file.path(plot_dir, "nonOutlier_Low_gene_spot_distribution.pdf"), width = 10, height = 6)
ggplot(nonOutlier_summ, aes(x = sample_id, y = n, fill = exclude_overlapping)) + geom_bar(stat="identity", position = position_dodge(preserve = "single")) +
facet_wrap(~precast_cluster, nrow = 3) + theme_classic() + ggtitle("Non-outlier (Low gene) spot distribution")
dev.off()

session_info()

