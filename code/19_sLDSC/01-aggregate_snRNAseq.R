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


res_dir <- here("processed-data", "19_sLDSC", "snRNA_seq", "input_files")
plot_dir <- here("plots", "19_sLDSC", "snRNA_seq")
dat_dir <- here("processed-data", "12_snRNA")
# Read in human snRNA-seq
sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
sce <- sce[ ,!sce$CellType.Final == "Neuron_Ambig"]

# Create pseudobulked profiles for each cell type
sce_pseudo <- aggregateAcrossCells(sce,DataFrame(cluster = colData(sce)$CellType.Final))

# Filter genes
# Based on expression
rowData(sce_pseudo)$high_expr_group_cluster <- filterByExpr(sce_pseudo, group = sce_pseudo$cluster)
sce_pseudo <- sce_pseudo[rowData(sce_pseudo)$high_expr_group_cluster, ]
sce_pseudo<-sce_pseudo[rowData(sce_pseudo)$gene_type=='protein_coding',]
sce_pseudo<-sce_pseudo[!duplicated(rowData(sce_pseudo)$gene_name),]
sce_pseudo <- sce_pseudo[!grepl("^MT-", rowData(sce_pseudo)$gene_name), ]
sce_pseudo <- sce_pseudo[!grepl("^RP[SL]", rowData(sce_pseudo)$gene_name), ]
rownames(sce_pseudo) <- rowData(sce_pseudo)$gene_name
# Normalize the expression
normd <- edgeR::cpm(edgeR::calcNormFactors(sce_pseudo))
colnames(normd)<-sce_pseudo$CellType.Final

write.table(normd, file.path(res_dir, 'snRNAseq_aggregated_cpm.tsv'), col.names = TRUE,
            row.names = TRUE, sep = "\t")

