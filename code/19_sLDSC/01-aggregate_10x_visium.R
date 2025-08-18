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

res_dir <- here("processed-data", "19_sLDSC", "10x_visium", "input_files")
raw_in_path <- here(
    "processed-data", "05_harmony_BayesSpace", "02-compute_QC_metrics", "spe_with_QC_metrics_hdf5"
)
spe <- loadHDF5SummarizedExperiment(raw_in_path)

# Add spatial domain information
clusters_resFile <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters", "precast_clusters.csv")
spe[["spatial_domains"]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(clusters_resFile), by = 'key') |>
    pull(cluster) |>
    as.factor()
# Remove spots with no PRECAST output
spe <- spe[ ,!is.na(spe[["spatial_domains"]])]
spe$spatial_domains <- as.character(spe$spatial_domains)
spe$spatial_domains <- gsub(" ", "_", spe$spatial_domains)
spe$spatial_domains <- gsub("/", "_", spe$spatial_domains)
spe$spatial_domains <- as.factor(spe$spatial_domains)

# Convert this to avoid issues with image data
counts <- assay(spe, "counts")
row_data <- rowData(spe)
col_data <- colData(spe)

spe2 <- SingleCellExperiment(
  assays = list(counts = as.matrix(counts)),
  rowData = row_data,
  colData = col_data
)

# Create pseudobulked profiles for each cell type
spe_pseudo <- aggregateAcrossCells(
    spe2, DataFrame(cluster = spe2[["spatial_domains"]]))

# Filter genes
# Based on expression
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$cluster)
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_cluster, ]
spe_pseudo<-spe_pseudo[rowData(spe_pseudo)$gene_type=='protein_coding',]
spe_pseudo<-spe_pseudo[!duplicated(rowData(spe_pseudo)$gene_name),]
spe_pseudo <- spe_pseudo[!grepl("^MT-", rowData(spe_pseudo)$gene_name), ]
spe_pseudo <- spe_pseudo[!grepl("^RP[SL]", rowData(spe_pseudo)$gene_name), ]
rownames(spe_pseudo) <- rowData(spe_pseudo)$gene_name

# Normalize the expression
normd <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo))
colnames(normd)<-spe_pseudo$spatial_domains

write.table(normd, file.path(res_dir, '10x_visium_aggregated_cpm.tsv'), col.names = TRUE,
            row.names = TRUE, sep = "\t")
