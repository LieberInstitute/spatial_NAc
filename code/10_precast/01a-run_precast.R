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
    c("nnSVG_type", "n", 1, "logical", "Use nnSVGs identified by controlling for PRECAST k = 2 clusters?"),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

k <- c(3:28)

spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5"
)

if (opt$nnSVG_type) {
    svg_path <- here(
        "processed-data", "05_harmony_BayesSpace", "07-run_nnSVG", "nnSVG_precast_out",
        "summary_across_samples.csv"
    )
    out_path <- here("processed-data", "10_precast", "nnSVG_precast")
    plot_dir <- here('plots', '10_precast', 'nnSVG_precast', 'BIC_select')
} else {
    svg_path <- here(
        "processed-data", "05_harmony_BayesSpace", "07-run_nnSVG", "nnSVG_out",
        "summary_across_samples.csv"
    )
        out_path <- here("processed-data", "10_precast", "nnSVG_default")
        plot_dir <- here('plots', '10_precast', 'nnSVG_default', 'BIC_select')
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
pre_obj <- AddParSetting(pre_obj, Sigma_equal = FALSE, verbose = TRUE, maxIter = 30)

#   Fit model
pre_obj <- PRECAST(pre_obj, K = k)
resList <- pre_obj@resList

saveRDS(pre_obj, file.path(out_path, "/PRECAST_BIC.rds"))
saveRDS(resList, file.path(out_path, "/resList_BIC.rds"))
session_info()
