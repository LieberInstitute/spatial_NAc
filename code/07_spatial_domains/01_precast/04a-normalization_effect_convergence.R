library(here)
library(PRECAST)
library(HDF5Array)
library(Seurat)
library(mclust)
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
library(ggpubr)

# Read and process the spe without genes filtered
spe_no_genes_filtered_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "02-compute_QC_metrics", "spe_with_QC_metrics_hdf5"
)
spe_no_genes_filtered <- loadHDF5SummarizedExperiment(spe_no_genes_filtered_dir)
spe_no_genes_filtered$low_umi <- spe_no_genes_filtered$sum_umi < 250
spe_no_genes_filtered$low_gene_edge_spot <- spe_no_genes_filtered$low_umi & (spe_no_genes_filtered$edge_distance < 6)
spe_no_genes_filtered <- spe_no_genes_filtered[ ,spe_no_genes_filtered$low_gene_edge_spot == FALSE]
spe_no_genes_filtered <- spe_no_genes_filtered[ ,spe_no_genes_filtered$local_outliers == "FALSE"]
# Filter genes with no expression in all spots
exprLogic <- counts(spe_no_genes_filtered) > 0
spe_no_genes_filtered$sample_id <- as.factor(spe_no_genes_filtered$sample_id)
nSpots_by_donor <- lapply(levels(spe_no_genes_filtered$sample_id), function(iSample){
    rowSums(exprLogic[ ,spe_no_genes_filtered$sample_id == iSample])
})
nSpots_by_donor <- do.call(cbind, nSpots_by_donor)
colnames(nSpots_by_donor) <- levels(spe_no_genes_filtered$sample_id)
select.genes <-  rownames(nSpots_by_donor)[rowSums(nSpots_by_donor == 0) == 0]
# Only select those genes which have non-zero expression in atleast 1 spot in each slide
spe_no_genes_filtered <- spe_no_genes_filtered[rownames(spe_no_genes_filtered) %in% select.genes, ]

# Read in the spe that has mito, ribo, and MALAT1 excluded
spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_dimRed_hdf5"
)
spe <- loadHDF5SummarizedExperiment(spe_dir)

# Read in SVGs
svg_path <- here(
     "processed-data", "05_harmony_BayesSpace", "07-run_nnSVG", "nnSVG_precast_out",
    "summary_across_samples.csv")

# Process the two spe's for PRECAST
#   PRECAST expects array coordinates in 'row' and 'col' columns
spe_no_genes_filtered$row <- spe_no_genes_filtered$array_row_transformed
spe_no_genes_filtered$col <- spe_no_genes_filtered$array_col_transformed
rownames(spe_no_genes_filtered) <- rowData(spe_no_genes_filtered)$gene_name
rownames(spe_no_genes_filtered) <- make.names(rownames(spe_no_genes_filtered), unique = TRUE)
seu_list_no_genes_filtered <- lapply(
    levels(spe_no_genes_filtered$donor),
    function(donor) {
        cat(donor, "\n")
        small_spe <- spe_no_genes_filtered[, spe_no_genes_filtered$donor == donor]
        cat(dim(small_spe), "\n")
        CreateSeuratObject(
            #   Bring into memory to greatly improve speed
            counts = as(assays(small_spe)$counts, "dgCMatrix"),
            meta.data = as.data.frame(colData(small_spe)),
            project = "spatialNAc"
        )
    }
)

spe$row <- spe$array_row_transformed
spe$col <- spe$array_col_transformed
rownames(spe) <- rowData(spe)$gene_name
rownames(spe) <- make.names(rownames(spe), unique = TRUE)
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

# Process SVGs
num_genes <- 2000
svgs_no_genes_filtered <- read.csv(svg_path) |>
    mutate(gene_name = rowData(spe_no_genes_filtered)[ match(gene_id, rowData(spe_no_genes_filtered)$gene_id), "gene_name"]) |>
    as_tibble() |>
    arrange(nnsvg_avg_rank_rank) |>
    slice_head(n = num_genes) |>
    pull(gene_name)
svgs <- read.csv(svg_path) |>
    mutate( gene_name = rowData(spe)[ match(gene_id, rowData(spe)$gene_id), "gene_name"]) |>
    as_tibble() |>
    arrange(nnsvg_avg_rank_rank) |>
    slice_head(n = num_genes) |>
    pull(gene_name)

# Create PRECAST objects
pre_obj_no_genes_filtered <- CreatePRECASTObject(
    seuList = seu_list_no_genes_filtered,
    selectGenesMethod = NULL,
    customGenelist = svgs_no_genes_filtered
    #   Using defaults for gene-filtering-related parameters. Though each donor
    #   consists of more spots than 1 typical Visium capture area (and would
    #   thus be expected to throw off the appropriateness of the defaults for
    #   'premin.spots', etc), we're using SVGs from nnSVG as input, and these
    #   genes already passed a similar reasonable expression cutoff:
    #   https://github.com/LieberInstitute/spatial_NAc/blob/61d1e198536a80bddca93017ea6eb8169af5d978/code/05_harmony_BayesSpace/05-run_nnSVG.R#L40-L45
)

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

pre_obj_no_genes_filtered <- AddAdjList(pre_obj_no_genes_filtered, platform = "Visium")
pre_obj <- AddAdjList(pre_obj, platform = "Visium")

# Fit the PRECAST model in both settings
s <- 1
pre_obj <- AddParSetting(
    pre_obj,
    Sigma_equal = FALSE, verbose = TRUE, maxIter = 30, coreNum_int=detectCores() - 1, seed = s)
pre_obj <- PRECAST(pre_obj, K = 6)

pre_obj_no_genes_filtered <- AddParSetting(pre_obj_no_genes_filtered, Sigma_equal = FALSE, verbose = TRUE, maxIter = 30, coreNum_int=detectCores() - 1, seed = s)
pre_obj_no_genes_filtered <- PRECAST(pre_obj_no_genes_filtered, K = 6)

resDir <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast")
saveRDS(pre_obj, file.path(resDir, "pre_obj_filtered_k6.rds"))
saveRDS(pre_obj_no_genes_filtered, file.path(resDir, "pre_obj_no_genes_filtered_k6.rds"))

ll_seq <- pre_obj@resList[[1]]$loglik_seq
ll_seq_no_genes_filtered <- pre_obj_no_genes_filtered@resList[[1]]$loglik_seq
ll_res <- data.frame("iteration" = c(1:30), "all_genes" = ll_seq_no_genes_filtered, "filtered_genes" = ll_seq)
ll_res <- reshape2::melt(ll_res, id.vars = "iteration")
colnames(ll_res) <- c("iteration", "gene_selection", "loglik")

plotDir <- here("plots", "07_spatial_domains", "01_precast", "nnSVG_precast")
pdf(file.path(plotDir, "comparing_convergence_gene_filtering.pdf"), height = 4, width = 5)
ggplot(ll_res[ll_res$iteration >= 2, ], aes(x = iteration, y = loglik, color = gene_selection)) + geom_point() + theme_classic() + ggtitle("K = 10")
dev.off()

# Compare to results without setting a limit on the number of iterations
pre_obj_final_path <- here::here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "random_start_1", "pre_obj_k10.rds")
pre_obj_final <- readRDS(pre_obj_final_path)

# Get results for intermediate and final preobj objects
pre_obj <- SelectModel(pre_obj, return_para_est=TRUE)
seuInt <- IntegrateSpaData(pre_obj, species = "Human")
#   Extract PRECAST results, clean up column names, and export to CSV
precast_results <- seuInt@meta.data |>
        rownames_to_column("key") |>
        as_tibble() |>
        select(-orig.ident) |>
        rename_with(~ sub("_PRE_CAST", "", .x))

pre_obj_final <- SelectModel(pre_obj_final, return_para_est=TRUE)
seuInt_final <- IntegrateSpaData(pre_obj_final, species = "Human")
#   Extract PRECAST results, clean up column names, and export to CSV
precast_results_final <- seuInt_final@meta.data |>
        rownames_to_column("key") |>
        as_tibble() |>
        select(-orig.ident) |>
        rename_with(~ sub("_PRE_CAST", "", .x))

precast_results <- precast_results |> select(c(key, cluster))
precast_results_final <- precast_results_final |> select(c(key, cluster))

temp <- colnames(spe)
colData(spe) <- colData(spe) |>
    as_tibble() |>
    left_join(precast_results, by = "key") |>
    DataFrame()
colnames(spe) <- temp

colnames(precast_results_final) <- c("key", "cluster_final")
temp <- colnames(spe)
colData(spe) <- colData(spe) |>
    as_tibble() |>
    left_join(precast_results_final, by = "key") |>
    DataFrame()
colnames(spe) <- temp

spot_plot(spe,sample_id = "Br2720", var_name = "cluster", is_discrete = TRUE, spatial = TRUE) + ggtitle("30 Iterations")+
guides(fill = guide_legend(override.aes = list(size = 5)))

spot_plot(spe,sample_id = "Br2720", var_name = "cluster_final", is_discrete = TRUE, spatial = TRUE) + ggtitle("Convergence")+
guides(fill = guide_legend(override.aes = list(size = 5)))