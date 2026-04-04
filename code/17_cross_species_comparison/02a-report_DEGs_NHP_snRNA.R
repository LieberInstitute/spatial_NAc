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

opt <- list()
opt$subset_neurons <- TRUE

if (is.null(opt$subset_neurons)) {
  opt$subset_neurons <- FALSE
}

print(opt$subset_neurons)

dat_dir <- here("processed-data", "12_snRNA")
ext_dir <- here("processed-data", "17_cross_species_comparison", "NHP_data")
reg_dir <- here("processed-data", "10_post_clustering_analysis", "02_spatial_registration_sn")

# Read human discovery dataset
sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
sce <- sce[, !sce$CellType.Final == "Neuron_Ambig"]

# Read NHP external dataset
if (opt$subset_neurons) {
  sce_ext <- readRDS(file = file.path(ext_dir, "GSE167920_Results_MSNs_processed_final.rds"))
} else {
  sce_ext <- readRDS(file = file.path(ext_dir, "GSE167920_Results_full_nuclei_processed_final.rds"))
}

if (opt$subset_neurons) {
  Idents(sce_ext) <- sce_ext$MSN_type
  Idents(sce_ext) <- gsub(Idents(sce_ext), pattern = "-", replacement = ".")
} else {
  Idents(sce_ext) <- sce_ext$cell_type
  Idents(sce_ext) <- gsub(Idents(sce_ext), pattern = "-", replacement = ".")
}

# Read human DEG results
sn_registration <- readRDS(file.path(reg_dir, "sn_cellType_registration.rds"))
sn_enrichment <- sn_registration[["enrichment"]]

if (opt$subset_neurons) {
  select_cell_types <- c(
    "DRD1_MSN_A", "DRD1_MSN_B", "DRD1_MSN_C", "DRD1_MSN_D",
    "DRD2_MSN_A", "DRD2_MSN_B"
  )
} else {
  select_cell_types <- c(
    "DRD1_MSN_A", "DRD1_MSN_B", "DRD1_MSN_C", "DRD1_MSN_D",
    "DRD2_MSN_A", "DRD2_MSN_B", "Excitatory", "Inh_A", "Inh_B",
    "Inh_C", "Inh_D", "Inh_E", "Inh_F", "Astrocyte_A", "Astrocyte_B",
    "Microglia", "OPC", "Oligo", "Ependymal", "Endothelial"
  )
}

sce <- sce[, sce$CellType.Final %in% select_cell_types]

# Build DEG table
DEG_Lists <- list()

for (i in select_cell_types) {
  df <- sn_enrichment[, grep(i, colnames(sn_enrichment)), drop = FALSE]
  colnames(df) <- gsub(paste0("_", i), "", colnames(df))
  df <- df[df$fdr < 0.05, , drop = FALSE]
  df <- df[df$logFC > 0, , drop = FALSE]
  df <- df[order(df$logFC, decreasing = TRUE), , drop = FALSE]
  df <- df[1:min(250, nrow(df)), , drop = FALSE]
  df$CellType <- i
  df$Gene_ID <- rownames(df)
  DEG_Lists[[i]] <- df
}

DEGs_df <- do.call(rbind, DEG_Lists)
rownames(DEGs_df) <- seq_len(nrow(DEGs_df))

cat("Rows in DEGs_df after combining cell-type DEG lists:", nrow(DEGs_df), "\n")
cat("Unique Gene_ID after combining cell-type DEG lists:", length(unique(DEGs_df$Gene_ID)), "\n")

# Add human gene symbols
geneData <- rowData(sce)
DEGs_df$Gene_symbol <- NA_character_

for (i in seq_len(nrow(DEGs_df))) {
  DEGs_df$Gene_symbol[i] <- geneData$gene_name[geneData$gene_id == DEGs_df$Gene_ID[i]][1]
}

cat("Unique human genes before removing missing symbols:", length(unique(DEGs_df$Gene_symbol)), "\n")

DEGs_df <- DEGs_df[!is.na(DEGs_df$Gene_symbol), , drop = FALSE]

cat("Rows after removing missing human symbols:", nrow(DEGs_df), "\n")
cat("Unique human genes after removing missing symbols:", length(unique(DEGs_df$Gene_symbol)), "\n")

# Add macaque ortholog symbols
Monkey_Orthos <- orthologs(
  genes = unique(DEGs_df$Gene_symbol),
  species = "Macaca mulatta"
)
colnames(Monkey_Orthos)[5] <- "Monkey_Symbol"

cat("Rows in Monkey_Orthos returned by babelgene:", nrow(Monkey_Orthos), "\n")
cat("Unique human symbols in Monkey_Orthos:", length(unique(Monkey_Orthos$human_symbol)), "\n")
cat("Unique macaque symbols in Monkey_Orthos:", length(unique(Monkey_Orthos$Monkey_Symbol)), "\n")

# Merge DEGs with ortholog results
DEGs_merged <- merge(
  x = DEGs_df,
  y = Monkey_Orthos,
  by.x = "Gene_symbol",
  by.y = "human_symbol"
)

cat("Rows after merging DEGs with macaque ortholog table:", nrow(DEGs_merged), "\n")
cat("Unique human genes after merge:", length(unique(DEGs_merged$Gene_symbol)), "\n")
cat("Unique macaque genes after merge:", length(unique(DEGs_merged$Monkey_Symbol)), "\n")

# Keep only mapped macaque symbols
DEGs_mapped <- DEGs_merged[!is.na(DEGs_merged$Monkey_Symbol), , drop = FALSE]

cat("Rows after removing missing macaque orthologs:", nrow(DEGs_mapped), "\n")
cat("Unique human genes after removing missing macaque orthologs:", length(unique(DEGs_mapped$Gene_symbol)), "\n")
cat("Unique macaque genes after removing missing macaque orthologs:", length(unique(DEGs_mapped$Monkey_Symbol)), "\n")

# True 1:1 ortholog filtering
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

# Require presence in macaque dataset
nhp_genes <- rownames(sce_ext[["RNA"]]$data)

DEGs_final <- DEGs_one2one[
  DEGs_one2one$Monkey_Symbol %in% nhp_genes,
  ,
  drop = FALSE
]

cat("Rows after requiring presence in macaque dataset:", nrow(DEGs_final), "\n")
cat("Unique human genes after requiring presence in macaque dataset:", length(unique(DEGs_final$Gene_symbol)), "\n")
cat("Unique macaque genes after requiring presence in macaque dataset:", length(unique(DEGs_final$Monkey_Symbol)), "\n")

# Summary of human gene counts
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

key_genes <- c(
  "DRD1", "DRD2", "TAC1", "ADORA2A", "PDYN",
  "PENK", "CALB1", "CPNE4", "RXFP1", "OPRM1"
)
check_df <- data.frame(
  Gene = key_genes,
  In_initial = key_genes %in% initial_human_genes,
  After_mapping = key_genes %in% mapped_human_genes,
  After_1to1 = key_genes %in% one2one_human_genes,
  In_final = key_genes %in% final_human_genes
)

print(check_df)