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

# Directories with results
res_dir <- here::here("processed-data", "23_functional_data_registration")
plot_dir <- here::here("plots", "23_functional_data_registration")

cocaine_acute <- read.csv(file.path(res_dir, "raw-data", "Cocaine_acute.csv"))
cocaine_chronic <- read.csv(file.path(res_dir, "raw-data", "Cocaine_repeated.csv"))
cocaine_all <- read.csv(file.path(res_dir, "raw-data", "Cocaine_all.csv"))

rat_human_ortholog <- read.delim("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/HOM_AllOrganism.rpt")
rat_orthologs <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "rat", ]
human_orthologs <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "human", ]
common_DB_class_keys <- intersect(rat_orthologs$DB.Class.Key, human_orthologs$DB.Class.Key)
rat_orthologs <- rat_orthologs[rat_orthologs$DB.Class.Key %in% common_DB_class_keys, ]
human_orthologs <- human_orthologs[human_orthologs$DB.Class.Key %in% common_DB_class_keys, ]
ortholog_list <- lapply(common_DB_class_keys, function(id){
rat_genes <- rat_orthologs$Symbol[rat_orthologs$DB.Class.Key == id]
human_genes <- human_orthologs$Symbol[human_orthologs$DB.Class.Key == id]
    list("rat_genes" = rat_genes, "human_genes" = human_genes)
})
names(ortholog_list) <- common_DB_class_keys
names(ortholog_list) <- common_DB_class_keys
length_rat_genes <- lapply(ortholog_list, function(iorg){
    length(iorg$rat_genes)
}) 
length_rat_genes <- unlist(length_rat_genes)

length_human_genes <- lapply(ortholog_list, function(iorg){
    length(iorg$human_genes)
}) 
length_human_genes <- unlist(length_human_genes)

ortholog_list <- ortholog_list[!((length_rat_genes == 0) | (length_human_genes == 0))]
length_rat_genes <- lapply(ortholog_list, function(iorg){
    length(iorg$rat_genes)
}) 
length_rat_genes <- unlist(length_rat_genes)

length_human_genes <- lapply(ortholog_list, function(iorg){
    length(iorg$human_genes)
}) 
length_human_genes <- unlist(length_human_genes)
ortholog_list_unique <- ortholog_list[(length_rat_genes == 1) & (length_human_genes == 1)]
rat_genes <- lapply(ortholog_list_unique, function(iorg){
    iorg$rat_genes
})
human_genes <- lapply(ortholog_list_unique, function(iorg){
    iorg$human_genes
})
rat_genes <- unlist(rat_genes)
human_genes <- unlist(human_genes)
orthologs_unique_df <- data.frame("JAX_id" = names(rat_genes), "rat_genes" = rat_genes, "human_genes" = human_genes)
rownames(orthologs_unique_df) <- c(1:dim(orthologs_unique_df)[1])

orthologs_df <- orthologs_unique_df[!duplicated(orthologs_unique_df$rat_genes), ]
orthologs_df <- orthologs_df[!duplicated(orthologs_df$human_genes), ]

# Now process the gene lists 
# Subset upregulated and downregulated genes
upregulated_acute <- subset(cocaine_acute, log2FoldChange > 0)
downregulated_acute <- subset(cocaine_acute, log2FoldChange < 0)
upregulated_chronic <- subset(cocaine_chronic, log2FoldChange > 0)
downregulated_chronic <- subset(cocaine_chronic, log2FoldChange < 0)
upregulated_all <- subset(cocaine_all, log2FoldChange > 0)
downregulated_all <- subset(cocaine_all, log2FoldChange < 0)

# Create lists of gene names grouped by Cell_type
upregulated_genes_acute <- split(upregulated_acute$Gene, upregulated_acute$Cell_type)
downregulated_genes_acute <- split(downregulated_acute$Gene, downregulated_acute$Cell_type)
upregulated_genes_chronic <- split(upregulated_chronic$Gene, upregulated_chronic$Cell_type)
downregulated_genes_chronic <- split(downregulated_chronic$Gene, downregulated_chronic$Cell_type)
upregulated_genes_all <- split(upregulated_all$Gene, upregulated_all$Cell_type)
downregulated_genes_all <- split(downregulated_all$Gene, downregulated_all$Cell_type)

# Convert the gene names from RAT to HUMAN using orthologs_df
convert_rat_to_human_genes <- function(deg_list, orthologs_df) {
  # Create a named vector for fast lookup: names are rat genes, values are human genes
  rat_to_human <- setNames(orthologs_df$human_genes, orthologs_df$rat_genes)
  
  # Apply the mapping to each list element
  human_deg_list <- lapply(deg_list, function(genes) {
    mapped <- rat_to_human[genes]
    # Remove NA (unmapped genes)
    mapped <- mapped[!is.na(mapped)]
    unname(mapped)
  })
  
  return(human_deg_list)
}
upregulated_genes_acute_human <- convert_rat_to_human_genes(upregulated_genes_acute, orthologs_df)
downregulated_genes_acute_human <- convert_rat_to_human_genes(downregulated_genes_acute, orthologs_df)
upregulated_genes_chronic_human <- convert_rat_to_human_genes(upregulated_genes_chronic, orthologs_df)
downregulated_genes_chronic_human <- convert_rat_to_human_genes(downregulated_genes_chronic, orthologs_df)
upregulated_genes_all_human <- convert_rat_to_human_genes(upregulated_genes_all, orthologs_df)
downregulated_genes_all_human <- convert_rat_to_human_genes(downregulated_genes_all, orthologs_df)
# Read in the spe data
raw_in_path <- here(
    "processed-data", "05_harmony_BayesSpace", "02-compute_QC_metrics", "spe_with_QC_metrics_hdf5"
)
spe <- loadHDF5SummarizedExperiment(raw_in_path)
rownames(spe) <- rowData(spe)$gene_name
gene_data <- rowData(spe)

# Convert gene names to gene IDs
convert_gene_names_to_ids <- function(gene_list, gene_data) {
  # Create mapping vector: gene_name -> gene_id
  name_to_id <- setNames(as.character(gene_data$gene_id), as.character(gene_data$gene_name))
  
  # Apply mapping
  lapply(gene_list, function(genes) {
    mapped <- name_to_id[genes]
    mapped <- mapped[!is.na(mapped)]  # remove unmapped
    unname(mapped)
  })
}
upregulated_genes_acute_ensembl <- convert_gene_names_to_ids(upregulated_genes_acute_human, gene_data)
downregulated_genes_acute_ensembl <- convert_gene_names_to_ids(downregulated_genes_acute_human, gene_data)
upregulated_genes_chronic_ensembl <- convert_gene_names_to_ids(upregulated_genes_chronic_human, gene_data)
downregulated_genes_chronic_ensembl <- convert_gene_names_to_ids(downregulated_genes_chronic_human, gene_data)
upregulated_genes_all_ensembl <- convert_gene_names_to_ids(upregulated_genes_all_human, gene_data)
downregulated_genes_all_ensembl <- convert_gene_names_to_ids(downregulated_genes_all_human, gene_data)

# Remove groups with fewer than 25 genes
upregulated_genes_acute_ensembl <- upregulated_genes_acute_ensembl[sapply(upregulated_genes_acute_ensembl, length) > 25]
downregulated_genes_acute_ensembl <- downregulated_genes_acute_ensembl[sapply(downregulated_genes_acute_ensembl, length) > 25]
upregulated_genes_chronic_ensembl <- upregulated_genes_chronic_ensembl[sapply(upregulated_genes_chronic_ensembl, length) > 25]
downregulated_genes_chronic_ensembl <- downregulated_genes_chronic_ensembl[sapply(downregulated_genes_chronic_ensembl, length) > 25]
upregulated_genes_all_ensembl <- upregulated_genes_all_ensembl[sapply(upregulated_genes_all_ensembl, length) > 25]
downregulated_genes_all_ensembl <- downregulated_genes_all_ensembl[sapply(downregulated_genes_all_ensembl, length) > 25]

# Read modeling results
modeling_results <- readRDS(file.path(res_dir, "modeling_results.rds"))

genes_all_ensembl <- list("Upregulated" = upregulated_genes_all_ensembl[["Drd1.MSN.1"]], 
                          "Downregulated" = downregulated_genes_all_ensembl[["Drd1.MSN.1"]])

cocaine_enrichment_all <- gene_set_enrichment(
    gene_list = genes_all_ensembl,
    fdr_cut = 0.05,
    modeling_results = modeling_results,
    model_type = "enrichment"
)

source(here("code", "23_functional_data_registration", "gene_set_enrichment_plot_complex.R"))
gene_count_col <- get_gene_list_count(genes_all_ensembl)
gene_count_row <- get_gene_enrichment_count(model_results = modeling_results, fdr_cut = 0.05)
cocaine_enrichment_all$test <- gsub("\\.", " ", cocaine_enrichment_all$test)
cocaine_enrichment_all$test[cocaine_enrichment_all$test == "Endothelial Ependymal"] <- "Endothelial/\nEpendymal"

rownames(gene_count_row) <- gsub("\\.", " ", rownames(gene_count_row))
rownames(gene_count_row)[rownames(gene_count_row) == "Endothelial Ependymal"] <- "Endothelial/\nEpendymal"
plot_dir <- here::here("plots", "23_functional_data_registration", "cocaine")
pdf(file.path(plot_dir, "Cocaine_all.pdf"), width = 5, height = 5)
gene_set_enrichment_plot_complex(
    enrichment = cocaine_enrichment_all,
    gene_count_col = gene_count_col,
    gene_count_row = gene_count_row,
    PThresh = 12,
    ORcut = 3,
    enrichOnly = FALSE,
    mypal = c(
        "white",
        grDevices::colorRampPalette(
            RColorBrewer::brewer.pal(9, "YlOrRd")
        )(50)
    ),
    anno_title_col = "Gene Set Size",
    anno_title_row = "Enriched Genes"
)
dev.off()