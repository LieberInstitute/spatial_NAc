library(SpatialExperiment)
library(SingleCellExperiment)
library(HDF5Array)
library(Seurat)
library(RColorBrewer)
library(spatialLIBD)
library(jaffelab)
library(here)
library(sessioninfo)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(fastTopics)
library(getopt)
library(edgeR)
library(scran)
library(scuttle)

# Read in the data
dat_dir <- here("code", "06_deploy_app", "spe_shiny")
spe <- loadHDF5SummarizedExperiment(dat_dir)

# Prepare expression data
expr_mat <- logcounts(spe)
geneData <- rowData(spe)
rownames(expr_mat) <- geneData$gene_id
# Only retain unique, protein coding autosomal genes
geneData <- geneData[!duplicated(geneData$gene_name), ]
geneData <- geneData[geneData$gene_type == "protein_coding", ]
geneData <- geneData[!grepl("^MT-", geneData$gene_name), ]
geneData <- geneData[!grepl("^RP[SL]", geneData$gene_name), ]

# Subset the expression matrix
expr_mat <- expr_mat[rownames(expr_mat) %in% geneData$gene_id, ]

plot_dir <- here("plots", "19_sLDSC", "meringue_factors")
res_dir <- here("processed-data", "19_sLDSC", "meringue_factors", "input_files")
dir.create(plot_dir, recursive = TRUE)
dir.create(res_dir, recursive = TRUE)

# Check MCP factors
plot_sample_order <- c("Br2743", "Br6432", "Br6423", "Br2720", "Br6471", "Br6522", "Br8492", "Br8325", "Br8667", "Br3942")
vis_grid_gene(spe, geneid = "meringue_cluster_1", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plot_dir, "meringue_cluster_1.pdf"))

vis_grid_gene(spe, geneid = "meringue_cluster_2", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plot_dir, "meringue_cluster_2.pdf"))

vis_grid_gene(spe, geneid = "meringue_cluster_3", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plot_dir, "meringue_cluster_3.pdf"))

vis_grid_gene(spe, geneid = "meringue_cluster_4", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plot_dir, "meringue_cluster_4.pdf"))

# Now compute the correlations
pattern_cols <- grep("^meringue_cluster_", colnames(colData(spe)), value = TRUE)

compute_cor_matrix <- function(expr, patterns) {
  cor_mat <- matrix(NA, nrow = nrow(expr), ncol = length(patterns))
  rownames(cor_mat) <- rownames(expr)
  colnames(cor_mat) <- patterns

  for (pat in patterns) {
    print(pat)
    pat_values <- colData(spe)[[pat]]
    non_na_idx <- which(!is.na(pat_values))
    expr_sub <- expr[, non_na_idx, drop = FALSE]
    pat_sub <- pat_values[non_na_idx]
    expr_sub <- as.matrix(expr_sub)
    cors <- apply(expr_sub, 1, function(x) cor(x, pat_sub, method = "pearson"))
    cor_mat[, pat] <- cors}
  cor_mat
}

cor_mat <- compute_cor_matrix(expr_mat, pattern_cols)
geneData <- geneData[match(rownames(cor_mat), geneData$gene_id), ]
rownames(cor_mat) <- geneData$gene_name

cor_mat <- cor_mat[rowSums(is.na(cor_mat)) == 0, ]

write.table(cor_mat, file.path(res_dir, 'meringue_correlations.tsv'), col.names = TRUE,
            row.names = TRUE, sep = "\t")