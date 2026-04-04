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
library(babelgene)
library(pheatmap)

spe_dir <- here("processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5")
subset_neurons <- FALSE
print(subset_neurons)

# Read in human 10x Visium
spe <- loadHDF5SummarizedExperiment(spe_dir)

# Read in NHP reference
ext_dir <- here("processed-data", "17_cross_species_comparison", "NHP_data")

if (subset_neurons) {
  sce_ext <- readRDS(file = file.path(ext_dir, "GSE167920_Results_MSNs_processed_final.rds"))
} else {
  sce_ext <- readRDS(file = file.path(ext_dir, "GSE167920_Results_full_nuclei_processed_final.rds"))
}

if (subset_neurons) {
  Idents(sce_ext) <- sce_ext$MSN_type
  Idents(sce_ext) <- gsub(Idents(sce_ext), pattern = "-", replacement = ".")
} else {
  Idents(sce_ext) <- sce_ext$cell_type
  Idents(sce_ext) <- gsub(Idents(sce_ext), pattern = "-", replacement = ".")
}

# Load human spatial pseudobulk markers
reg_dir <- here(
  "processed-data", "10_post_clustering_analysis", "01_pseudobulk_markers",
  "01_precast", "pseudobulk_capture_area", "final_clusters"
)

load(file.path(reg_dir, "model_results_precast_clusters.Rdata"))
sp_enrichment <- modeling_results$enrichment
rm(modeling_results)

if (subset_neurons) {
  select_domains <- c("MSN.1", "MSN.2", "MSN.3", "D1.islands")
} else {
  select_domains <- c("MSN.1", "MSN.2", "MSN.3", "D1.islands", "Excitatory", "Inhibitory", "WM", "Endothelial.Ependymal")
}

# Add PRECAST labels
clusters_file <- here(
  "processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast",
  "final_clusters", "precast_clusters.csv"
)

spe[["precast_clusters"]] <- colData(spe) |>
  as_tibble() |>
  left_join(read.csv(clusters_file), by = "key") |>
  pull(cluster) |>
  as.factor()

spe <- spe[, !is.na(spe[["precast_clusters"]])]
spe[["precast_clusters"]] <- gsub(" ", ".", spe[["precast_clusters"]])
spe[["precast_clusters"]] <- gsub("/", ".", spe[["precast_clusters"]])
spe <- spe[, spe$precast_clusters %in% select_domains]

# Build DEG table
DEG_Lists <- list()

for (i in select_domains) {
  df <- sp_enrichment[, grep(i, colnames(sp_enrichment)), drop = FALSE]
  colnames(df) <- gsub(paste0("_", i), "", colnames(df))
  df <- df[df$fdr < 0.05, , drop = FALSE]
  df <- df[df$logFC > 0, , drop = FALSE]
  df <- df[order(df$logFC, decreasing = TRUE), , drop = FALSE]
  df <- df[1:min(250, dim(df)[1]), , drop = FALSE]
  df$CellType <- i
  df$Gene_ID <- rownames(df)
  DEG_Lists[[i]] <- df
}

DEGs_df <- do.call(rbind, DEG_Lists)
rownames(DEGs_df) <- seq_len(nrow(DEGs_df))

cat("Rows in DEGs_df after combining domain DEG lists:", nrow(DEGs_df), "\n")
cat("Unique Gene_ID after combining domain DEG lists:", length(unique(DEGs_df$Gene_ID)), "\n")

# Add human symbols
geneData <- rowData(spe)
DEGs_df$Gene_symbol <- NA_character_
for (i in seq_len(nrow(DEGs_df))) {
  DEGs_df$Gene_symbol[i] <- geneData$gene_name[geneData$gene_id == DEGs_df$Gene_ID[i]][1]
}

cat("Unique human genes before removing missing symbols:", length(unique(DEGs_df$Gene_symbol)), "\n")

DEGs_df <- DEGs_df[!is.na(DEGs_df$Gene_symbol), , drop = FALSE]

cat("Rows after removing missing human symbols:", nrow(DEGs_df), "\n")
cat("Unique human genes after removing missing symbols:", length(unique(DEGs_df$Gene_symbol)), "\n")

# Map human -> macaque orthologs
Monkey_Orthos <- orthologs(
  genes = unique(DEGs_df$Gene_symbol),
  species = "Macaca mulatta"
)
colnames(Monkey_Orthos)[5] <- "Monkey_Symbol"

cat("Rows in Monkey_Orthos returned by babelgene:", nrow(Monkey_Orthos), "\n")
cat("Unique human symbols in Monkey_Orthos:", length(unique(Monkey_Orthos$human_symbol)), "\n")
cat("Unique macaque symbols in Monkey_Orthos:", length(unique(Monkey_Orthos$Monkey_Symbol)), "\n")

DEGs_merged <- merge(
  x = DEGs_df,
  y = Monkey_Orthos,
  by.x = "Gene_symbol",
  by.y = "human_symbol"
)

cat("Rows after merging DEGs with macaque ortholog table:", nrow(DEGs_merged), "\n")
cat("Unique human genes after merge:", length(unique(DEGs_merged$Gene_symbol)), "\n")
cat("Unique macaque genes after merge:", length(unique(DEGs_merged$Monkey_Symbol)), "\n")

DEGs_mapped <- DEGs_merged[!is.na(DEGs_merged$Monkey_Symbol), , drop = FALSE]

cat("Rows after removing missing macaque orthologs:", nrow(DEGs_mapped), "\n")
cat("Unique human genes after removing missing macaque orthologs:", length(unique(DEGs_mapped$Gene_symbol)), "\n")
cat("Unique macaque genes after removing missing macaque orthologs:", length(unique(DEGs_mapped$Monkey_Symbol)), "\n")

# True 1:1 filtering
ortholog_pairs <- unique(DEGs_mapped[, c("Gene_symbol", "Monkey_Symbol")])

cat("Unique human-macaque ortholog pairs:", nrow(ortholog_pairs), "\n")

human_counts <- table(ortholog_pairs$Gene_symbol)
nhp_counts <- table(ortholog_pairs$Monkey_Symbol)

valid_pairs <- ortholog_pairs[
  ortholog_pairs$Gene_symbol %in% names(human_counts[human_counts == 1]) &
    ortholog_pairs$Monkey_Symbol %in% names(nhp_counts[nhp_counts == 1]),
  ,
  drop = FALSE
]

cat("Unique human genes in 1:1 ortholog pairs:", length(unique(valid_pairs$Gene_symbol)), "\n")
cat("Unique macaque genes in 1:1 ortholog pairs:", length(unique(valid_pairs$Monkey_Symbol)), "\n")
cat("Unique 1:1 ortholog pairs retained:", nrow(valid_pairs), "\n")

DEGs_one2one <- DEGs_mapped[
  DEGs_mapped$Gene_symbol %in% valid_pairs$Gene_symbol &
    DEGs_mapped$Monkey_Symbol %in% valid_pairs$Monkey_Symbol,
  ,
  drop = FALSE
]

cat("Rows after applying 1:1 ortholog filtering:", nrow(DEGs_one2one), "\n")
cat("Unique human genes after applying 1:1 ortholog filtering:", length(unique(DEGs_one2one$Gene_symbol)), "\n")
cat("Unique macaque genes after applying 1:1 ortholog filtering:", length(unique(DEGs_one2one$Monkey_Symbol)), "\n")

# Require presence in NHP dataset
nhp_genes <- rownames(sce_ext[["RNA"]]$data)

DEGs_final <- DEGs_one2one[
  DEGs_one2one$Monkey_Symbol %in% nhp_genes,
  ,
  drop = FALSE
]

cat("Rows after requiring presence in macaque dataset:", nrow(DEGs_final), "\n")
cat("Unique human genes after requiring presence in macaque dataset:", length(unique(DEGs_final$Gene_symbol)), "\n")
cat("Unique macaque genes after requiring presence in macaque dataset:", length(unique(DEGs_final$Monkey_Symbol)), "\n")

# Summary for table
initial_human_genes <- unique(DEGs_df$Gene_symbol)
mapped_human_genes <- unique(DEGs_mapped$Gene_symbol)
one2one_human_genes <- unique(DEGs_one2one$Gene_symbol)
final_human_genes <- unique(DEGs_final$Gene_symbol)

cat("\n--- Summary of human gene counts ---\n")
cat("Initial unique human genes:", length(initial_human_genes), "\n")
cat("After mapped macaque ortholog filter:", length(mapped_human_genes), "\n")
cat("After true 1:1 ortholog filter:", length(one2one_human_genes), "\n")
cat("After requiring macaque dataset expression:", length(final_human_genes), "\n")

cat("\n--- Genes lost at each step ---\n")
cat("Lost when requiring mapped macaque ortholog:", length(setdiff(initial_human_genes, mapped_human_genes)), "\n")
cat("Lost when requiring true 1:1 ortholog:", length(setdiff(mapped_human_genes, one2one_human_genes)), "\n")
cat("Lost when requiring macaque dataset expression:", length(setdiff(one2one_human_genes, final_human_genes)), "\n")

# Key-marker check
key_genes <- c("DRD1", "DRD2", "TAC1", "ADORA2A", "PDYN", "PENK", "CALB1", "CPNE4", "RXFP1", "OPRM1")

check_df <- data.frame(
  Gene = key_genes,
  In_initial = key_genes %in% initial_human_genes,
  After_mapping = key_genes %in% mapped_human_genes,
  After_1to1 = key_genes %in% one2one_human_genes,
  In_final = key_genes %in% final_human_genes
)

print(check_df)