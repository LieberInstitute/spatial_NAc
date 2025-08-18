rm(list = ls())
library(Seurat)
library(Matrix)
library(CoGAPS)
library(here)
library(ggplot2)
library(pheatmap)
library(HDF5Array)
library(SpatialExperiment)
library(spatialLIBD)
library(scran)
library(scuttle)
library(SingleCellExperiment)
library(RcppML)
library(tidyverse)
library(edgeR)
library(scater)
library(scran)
library(dplyr)
library(gridExtra)
library(ggforce)

# Create a directory to store the results
res_dir <- here::here("processed-data", "23_functional_data_registration")
plot_dir <- here::here("plots", "23_functional_data_registration")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Read in the SPE data
raw_in_path <- here(
    "processed-data", "05_harmony_BayesSpace", "02-compute_QC_metrics", "spe_with_QC_metrics_hdf5"
)
spe <- loadHDF5SummarizedExperiment(raw_in_path)
rownames(spe) <- rowData(spe)$gene_name
rownames(spe) <- make.names(rownames(spe), unique = TRUE)
spe$low_umi <- spe$sum_umi < 250
spe$low_gene <- spe$sum_gene < 250
spe$low_gene_edge_spot <- spe$low_umi & spe$edge_distance < 6
spe <- spe[ ,spe$low_gene_edge_spot == "FALSE"]
spe <- spe[ ,spe$local_outliers == "FALSE"]
exprLogic <- counts(spe) > 0
spe$sample_id <- as.factor(spe$sample_id)
nSpots_by_donor <- lapply(levels(spe$sample_id), function(iSample){
    rowSums(exprLogic[ ,spe$sample_id == iSample])
})
nSpots_by_donor <- do.call(cbind, nSpots_by_donor)
colnames(nSpots_by_donor) <- levels(spe$sample_id)
select.genes <-  rownames(nSpots_by_donor)[rowSums(nSpots_by_donor == 0) == 0]
# Only select those genes which have non-zero expression in atleast 1 spot in each slide
spe <- spe[rownames(spe) %in% select.genes, ]
spe <- logNormCounts(spe)

# Compute spatial domain registration
# First pseudobulk
cluster_dir <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters")
cluster_file <- "precast_clusters.csv"
agg_level <- "sample_id_original"
clusters_resFile <- file.path(cluster_dir, cluster_file)

spe[["spatial_domains"]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(clusters_resFile), by = 'key') |>
    pull(cluster) |>
    as.factor()
spe <- spe[ ,!is.na(spe[["spatial_domains"]])]
counts(spe) <- as(realize(counts(spe)), "dgCMatrix")

# Computing pseudobulk for registration
spe_pseudo <- aggregateAcrossCells(
    spe, DataFrame(cluster = spe[["spatial_domains"]], sample_id = colData(spe)$sample_id_original))
colnames(spe_pseudo) <- paste0(spe_pseudo$sample_id_original, "_", spe_pseudo[["spatial_domains"]])

colData(spe_pseudo) = colData(spe_pseudo)[
    , sort(c("sample_id", "sample_id_original", "slide_num", "in_tissue", "slide_num", "donor",  "Age", "Sex", "spatial_domains", "Diagnosis", "ncells"))]
colData(spe_pseudo) <- colData(spe_pseudo)[ , unique(colnames(colData(spe_pseudo)))]

spe_pseudo$sample_id_original <- as.character(spe_pseudo$sample_id_original)
spe_pseudo$sample_id_original <- make.names(spe_pseudo$sample_id_original)
spe_pseudo$sample_id_original <- factor(spe_pseudo$sample_id_original)

spe_pseudo$slide_num <- as.character(spe_pseudo$slide_num)
spe_pseudo$slide_num <- make.names(spe_pseudo$slide_num)
spe_pseudo$slide_num <- factor(spe_pseudo$slide_num)

spe_pseudo[["spatial_domains"]] <- as.character(spe_pseudo[["spatial_domains"]])
spe_pseudo[["spatial_domains"]] <- make.names(spe_pseudo[["spatial_domains"]])
spe_pseudo$spatial_domains <- factor(spe_pseudo$spatial_domains)

# Filter samples with too few nCells
print(dim(spe_pseudo))
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 50]
print(dim(spe_pseudo))

rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id_original)
rowData(spe_pseudo)$high_expr_group_domain <- filterByExpr(spe_pseudo, group = colData(spe_pseudo)[["spatial_domains"]])
summary(rowData(spe_pseudo)$high_expr_group_sample_id)
summary(rowData(spe_pseudo)$high_expr_group_domain)
with(rowData(spe_pseudo), table(high_expr_group_sample_id, high_expr_group_domain))
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_domain, ]
dim(spe_pseudo)

# Process the spe count data
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)
stopifnot(identical(rownames(x), rownames(spe_pseudo)))
dimnames(x) <- dimnames(spe_pseudo)
logcounts(spe_pseudo) <- x
spe_pseudo <- scuttle::addPerCellQC(spe_pseudo)
spe_pseudo$det_out <- spe_pseudo$detected < 6000
pdf(file = file.path(plot_dir, "detected_genes_domain.pdf"), width = 10, height = 6)
par(mar = c(10, 5, 2, 2))
boxplot(detected ~ colData(spe_pseudo)[["spatial_domains"]], data = colData(spe_pseudo), xlab = "", ylab = "Number of detected genes", las =2)
dev.off()
spe_pseudo<-spe_pseudo[,!spe_pseudo$det_out]
dim(spe_pseudo)
rm(x)

outlier_clusters <- names(table(spe_pseudo[["spatial_domains"]]))[table(spe_pseudo[["spatial_domains"]]) < 10] 
spe_pseudo <- spe_pseudo[ ,!as.character(spe_pseudo[["spatial_domains"]]) %in% outlier_clusters]
spe_pseudo[["spatial_domains"]] <- factor(as.character(spe_pseudo[["spatial_domains"]]))

# Print final spe dimension
print("Final SPE dimensions")
dim(spe_pseudo)

#run PCA
set.seed(12141)
#runPCA(spe_pseudo)
pca <- prcomp(t(assays(spe_pseudo)$logcounts))

message(Sys.time(), " % of variance explained by the PCs:")
metadata(spe_pseudo)
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca))
metadata(spe_pseudo)

pca_pseudo <- pca$x
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)
pdf(file = file.path(plot_dir, "pseudobulk_PCA_visium_final.pdf"), width = 8, height = 8)
plotPCA(spe_pseudo, colour_by = "sum", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "detected", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "ncells", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sample_id_original", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "spatial_domains", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
dev.off()

# Run DE model
# Rescale the feature age
spe_pseudo$age_scaled <- scales::rescale(spe_pseudo$Age,to=c(0,1))
if(!is.factor(spe_pseudo$Sex)){
    spe_pseudo$Sex <- factor(spe_pseudo$Sex, levels = c("F", "M"))
}
if(!is.factor(spe_pseudo$slide_num)){
    spe_pseudo$slide_num <- factor(spe_pseudo$slide_num, levels = unique(spe_pseudo$slide_num))
}
if(!is.factor(spe_pseudo$donor)){
    spe_pseudo$donor <- factor(spe_pseudo$donor, levels = unique(spe_pseudo$donor))
}
covars <- c('Sex', 'slide_num')
mod<-registration_model(spe_pseudo, covars = covars, var_registration = "spatial_domains")
cors<-registration_block_cor(spe_pseudo, registration_model = mod,  var_sample_id = "sample_id_original")
ensembl_col = 'gene_id'
symbol_col = 'gene_name'
results_enrichment <-registration_stats_enrichment(
    spe_pseudo,
    block_cor=cors,
    covars = covars,
    var_registration = "spatial_domains",
    var_sample_id = "sample_id_original",
    gene_ensembl = ensembl_col,
    gene_name = symbol_col
)

results_anova <-registration_stats_anova(
    spe_pseudo,
    block_cor=cors,
    covars = covars,
    var_registration = "spatial_domains",
    var_sample_id = "sample_id_original",
    gene_ensembl = ensembl_col,
    gene_name = symbol_col, 
    suffix = "all"
)

results_pairwise <-registration_stats_pairwise(
    spe_pseudo,
    block_cor=cors,
    registration_model=mod,
    var_registration = "spatial_domains",
    var_sample_id = "sample_id_original",
    gene_ensembl = ensembl_col,
    gene_name = symbol_col
)

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_enrichment,
    "pairwise" = results_pairwise
)

saveRDS(modeling_results, file.path(res_dir, "modeling_results.rds"))