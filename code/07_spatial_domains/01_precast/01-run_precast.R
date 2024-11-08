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
library(ggsci)
library(ggpubr)
library(getopt)
library(scran)
library(scater)
library(spatialNAcUtils)
library(spatialLIBD)

spec <- matrix(
    c(
        "nnSVG_type", "n", 1, "logical", "Use nnSVGs identified by controlling for PRECAST k = 2 clusters?",
        "specify_k", "k", 1, "logical", "Pre-specify number of clusters?", 
        "random_seed", "r", 1, "integer", "Random start number"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)
random_start_test <- TRUE

if(opt$specify_k){
    k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}else{
    k <- c(3:28)
}

spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_dimRed_hdf5"
)

if (opt$nnSVG_type) {
    svg_path <- here(
        "processed-data", "05_harmony_BayesSpace", "07-run_nnSVG", "nnSVG_precast_out",
        "summary_across_samples.csv"
    )
    if(opt$specify_k){
        if(random_start_test){
            out_path <- here("processed-data","07_spatial_domains", "01_precast", "nnSVG_precast", paste0("random_start_", opt$random_seed))
            plot_dir <- here("plots", "07_spatial_domains", "01_precast", "nnSVG_precast", paste0("random_start_", opt$random_seed), paste0("cluster_k_", k))
        }else{
        out_path <- here("processed-data","07_spatial_domains", "01_precast", "nnSVG_precast")
             plot_dir <- here('plots',"07_spatial_domains", '01_precast', 'nnSVG_precast', paste0("cluster_k_", k))
        }
       }
    else{
        out_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast")
        plot_dir <- here('plots', '07_spatial_domains', '01_precast', 'nnSVG_precast', 'BIC_select')
    }
} else {
    svg_path <- here(
        "processed-data", "05_harmony_BayesSpace", "07-run_nnSVG", "nnSVG_out",
        "summary_across_samples.csv"
    )
    if(opt$specify_k){
        if(random_start_test){
            out_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_default", paste0("random_start_", opt$random_seed))
            plot_dir <- here("plots", "07_spatial_domains", "01_precast", "nnSVG_precast", paste0("random_start_", opt$random_seed), paste0("cluster_k_", k))
        }else{
        out_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_default")
             plot_dir <- here('plots', '07_spatial_domains', '01_precast', 'nnSVG_precast', paste0("cluster_k_", k))
        }} 
    else {
        out_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_default")
        plot_dir <- here('plots', '07_spatial_domains', '01_precast', 'nnSVG_default', 'BIC_select')
    }
}
print("Output saved at:")
print(out_path)
print("Spatially variable genes obtained from:")
print(svg_path)
print("Plots saved at:")
print(plot_dir)
num_genes <- 2000

set.seed(1)
dir.create(dirname(out_path), showWarnings = FALSE)
dir.create(dirname(plot_dir), showWarnings = FALSE)

spe <- loadHDF5SummarizedExperiment(spe_dir)

#   PRECAST expects array coordinates in 'row' and 'col' columns
spe$row <- spe$array_row_transformed
spe$col <- spe$array_col_transformed

rownames(spe) <- rowData(spe)$gene_name
rownames(spe) <- make.names(rownames(spe), unique = TRUE)

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

svgs <- read.csv(svg_path) |>
    mutate( gene_name = rowData(spe)[ match(gene_id, rowData(spe)$gene_id), "gene_name"]) |>
    as_tibble() |>
    arrange(nnsvg_avg_rank_rank) |>
    slice_head(n = num_genes) |>
    pull(gene_name)

pre_obj <- CreatePRECASTObject(
    seuList = seu_list,
    selectGenesMethod = NULL,
    customGenelist = svgs
    #   Using defaults for gene-filtering-related parameters. Though each donor
    #   consists of more spots than 1 typical Visium capture area (and would
    #   thus be expected to throw off the appropriateness of the defaults for
    #   'premin.spots', etc), we're using SVGs from nnSVG as input, and these
    #   genes already passed a similar reasonable expression cutoff:
    #   https://github.com/LieberInstitute/spatial_NAc/blob/61d1e198536a80bddca93017ea6eb8169af5d978/code/05_harmony_BayesSpace/05-run_nnSVG.R#L40-L45
)

#   Setting platform to "Visium" just means to use array indices, which should
#   work fine despite the abnormal/ "artificial" capture area we've created by
#   stitching
pre_obj <- AddAdjList(pre_obj, platform = "Visium")

#   Following https://feiyoung.github.io/PRECAST/articles/PRECAST.BreastCancer.html,
#   which involves overriding some default values, though the implications are not
#   documented
if(random_start_test){
    if(opt$random_seed == 1){
        s <- 1
    }else{
        if(opt$random_seed == 2){
            s <- 1234
        }else{
            if(opt$random_seed == 3){
                s <- 1234567
            }else{
                if(opt$random_seed == 4){
                    s <- 100000001
                }else{
                    if(opt$random_seed == 5){
                        s <- 98989898
                    }else{
                        stop("Invalid random start provided")
                    }
                }
            }
        }
    }
    
    cat("Random seed = ", "\n")
    cat(s, "\n")
    cat("Out path = ", "\n")
    cat(paste0(out_path, "/pre_obj_k", k, ".rds"), "\n")
    pre_obj <- AddParSetting(
    pre_obj,
    Sigma_equal = FALSE, verbose = TRUE, maxIter = 500, coreNum_int=detectCores() - 1, seed = s)
    pre_obj <- PRECAST(pre_obj, K = k)
    # Save the output
    saveRDS(pre_obj, paste0(out_path, "/pre_obj_k", k, ".rds"))
    pre_obj <- SelectModel(pre_obj, return_para_est=TRUE)
    seuInt <- IntegrateSpaData(pre_obj, species = "Human")
    #   Extract PRECAST results, clean up column names, and export to CSV
    precast_results <- seuInt@meta.data |>
        rownames_to_column("key") |>
        as_tibble() |>
        select(-orig.ident) |>
        rename_with(~ sub("_PRE_CAST", "", .x))

    write_csv(precast_results, file.path(out_path, sprintf("PRECAST_k%s.csv", k)))

    cols_cluster <- get_palette(palette = "default", k)

    seuInt <- RunTSNE(seuInt, reduction = "PRECAST", tSNE.method = "FIt-SNE")

    p1 <- dimPlot(seuInt, item = "cluster", point_size = 0.5, font_family = "serif", cols = cols_cluster,
    border_col = "gray10", nrow.legend = 14, legend_pos = "right")

    cols_batch <- chooseColors(palettes_name = "Classic 20", n_colors = 10, plot_colors = TRUE)
    p2 <- dimPlot(seuInt, item = "batch", point_size = 0.5, font_family = "serif", cols = cols_batch,
    border_col = "gray10", nrow.legend = 14, legend_pos = "right")

    pdf(file.path(plot_dir, "/precast_tSNE.pdf"), width = 12, height = 8)
    plot_grid(p1, p2, ncol = 2)
    dev.off()

    dat_deg <- FindAllMarkers(seuInt, logfc.threshold = 0.1, only.pos = TRUE)
    write_csv(dat_deg, file.path(out_path, sprintf("PRECAST_k%s_marker_genes.csv", k)))
    n <- 10
    dat_deg %>%
        group_by(cluster) %>%
        top_n(n = n, wt = avg_log2FC) -> top_markers

    pdf(file.path(plot_dir, "/precast_cluster_markers.pdf"), width = 15, height = 7)
    DotPlot(seuInt, features = unique(top_markers$gene), col.min = 0, col.max = 1) + theme(axis.text.x = element_text(angle = 45,
        hjust = 1, size = 8))
    dev.off()

    precast_cluster <- colData(spe) |>
            as_tibble() |>
            left_join(precast_results, by = "key") |>
            pull(cluster)

    spe$precast_cluster <- precast_cluster
    # Extract PRECAST results 
    pdf(file.path(plot_dir, "precast_clusters.pdf"), width = 8, height = 8)
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

    pdf(file.path(plot_dir, "GLMPCA_precast_clusters.pdf"), width = 6, height = 6)
    plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 2, colour_by = "precast_cluster")
    dev.off()

    pdf(file.path(plot_dir, "PCA_precast_clusters.pdf"), width = 6, height = 6)
    plotReducedDim(spe, dimred = "PCA_p2", ncomponents = 2, colour_by = "precast_cluster")
    dev.off()

    saveRDS(seuInt, file.path(out_path, sprintf("/PRECAST_k%s_integrated.rds", k)))
}else{
    s <- 1
    pre_obj <- AddParSetting(
    pre_obj,
    Sigma_equal = FALSE, verbose = TRUE, maxIter = 500, coreNum_int=detectCores() - 1, seed = s)
    pre_obj <- PRECAST(pre_obj, K = k)
    resList <- pre_obj@resList
    pre_obj <- SelectModel(pre_obj, return_para_est=TRUE)
    pre_obj <- IntegrateSpaData(pre_obj, species = "Human")
    #   Extract PRECAST results, clean up column names, and export to CSV
    precast_results <- pre_obj@meta.data |>
        rownames_to_column("key") |>
        as_tibble() |>
        select(-orig.ident) |>
        rename_with(~ sub("_PRE_CAST", "", .x))

    write_csv(precast_results, file.path(out_path, sprintf("PRECAST_k%s.csv", k)))

    cols_cluster <- get_palette(palette = "default", k)

    pre_obj <- RunTSNE(pre_obj, reduction = "PRECAST", tSNE.method = "FIt-SNE")

    p1 <- dimPlot(pre_obj, item = "cluster", point_size = 0.5, font_family = "serif", cols = cols_cluster,
    border_col = "gray10", nrow.legend = 14, legend_pos = "right")

    cols_batch <- chooseColors(palettes_name = "Classic 20", n_colors = 10, plot_colors = TRUE)
    p2 <- dimPlot(pre_obj, item = "batch", point_size = 0.5, font_family = "serif", cols = cols_batch,
    border_col = "gray10", nrow.legend = 14, legend_pos = "right")

    pdf(file.path(plot_dir, "/precast_tSNE.pdf"), width = 12, height = 8)
    plot_grid(p1, p2, ncol = 2)
    dev.off()

    dat_deg <- FindAllMarkers(pre_obj, logfc.threshold = 0.1, only.pos = TRUE)
    write_csv(dat_deg, file.path(out_path, sprintf("PRECAST_k%s_marker_genes.csv", k)))
    n <- 10
    dat_deg %>%
        group_by(cluster) %>%
        top_n(n = n, wt = avg_log2FC) -> top_markers

    pdf(file.path(plot_dir, "/precast_cluster_markers.pdf"), width = 15, height = 7)
    DotPlot(pre_obj, features = unique(top_markers$gene), col.min = 0, col.max = 1) + theme(axis.text.x = element_text(angle = 45,
        hjust = 1, size = 8))
    dev.off()

    precast_cluster <- colData(spe) |>
            as_tibble() |>
            left_join(precast_results, by = "key") |>
            pull(cluster)

    spe$precast_cluster <- precast_cluster
    # Extract PRECAST results 
    pdf(file.path(plot_dir, "precast_clusters.pdf"), width = 8, height = 8)
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

    pdf(file.path(plot_dir, "GLMPCA_precast_clusters.pdf"), width = 6, height = 6)
    plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 2, colour_by = "precast_cluster")
    dev.off()

    pdf(file.path(plot_dir, "PCA_precast_clusters.pdf"), width = 6, height = 6)
    plotReducedDim(spe, dimred = "PCA_p2", ncomponents = 2, colour_by = "precast_cluster")
    dev.off()

    saveRDS(pre_obj, file.path(out_path, sprintf("/PRECAST_k%s_integrated.rds", k)))
}
 session_info()






