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
subset_neurons <- TRUE

print(subset_neurons)

# Read in human 10x Visium
spe <- loadHDF5SummarizedExperiment(spe_dir)

# Read in rodent reference
dat_dir <- here("processed-data", "12_snRNA")
sce_ext <- readRDS(file = file.path(dat_dir, "NAc_Combo_Integrated.RDS"))
sce_ext_neurons <- subset(
  sce_ext,
  idents = c(
    "Drd1-MSN-1", "Drd1-MSN-2", "Drd2-MSN-1", "Drd2-MSN-2",
    "Drd3-MSN", "Grm8-MSN", "GABAergic", "Chat-Interneuron",
    "Pvalb-Interneuron", "Sst-Interneuron", "Glutamatergic"
  )
)

Idents(sce_ext_neurons) <- gsub(Idents(sce_ext_neurons), pattern = "-", replacement = ".")
Idents(sce_ext) <- gsub(Idents(sce_ext), pattern = "-", replacement = ".")

# Read spatial marker results
reg_dir <- here(
  "processed-data", "10_post_clustering_analysis", "01_pseudobulk_markers",
  "01_precast", "pseudobulk_capture_area", "final_clusters"
)
load(file.path(reg_dir, "model_results_precast_clusters.Rdata"))
sp_enrichment <- modeling_results$enrichment
rm(modeling_results)

if (subset_neurons) {
  select_domains <- c("MSN.1", "MSN.2", "MSN.3", "D1.islands", "Excitatory", "Inhibitory")
  rat_genes <- rownames(sce_ext_neurons)
} else {
  select_domains <- c("MSN.1", "MSN.2", "MSN.3", "D1.islands", "Excitatory", "Inhibitory", "WM", "Endothelial.Ependymal")
  rat_genes <- rownames(sce_ext)
}

# Add PRECAST results
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

# Add rat orthologs
rat_human_ortholog <- read.delim("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/HOM_AllOrganism.rpt")
hom_rat <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "rat", ]
hom_hs <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "human", ]

DEGs_df$Rat_Gene_Symbol <- NA_character_

for (i in seq_len(nrow(DEGs_df))) {
  human_keys <- hom_hs[hom_hs$Symbol %in% DEGs_df[i, "Gene_symbol"], "DB.Class.Key"]
  rat_symbols <- unique(hom_rat[hom_rat$DB.Class.Key %in% human_keys, "Symbol"])

  if (length(rat_symbols) == 1) {
    DEGs_df[i, "Rat_Gene_Symbol"] <- rat_symbols
  }
}

cat("Unique human genes before removing missing rat orthologs:", length(unique(DEGs_df$Gene_symbol)), "\n")
cat("Unique rat ortholog symbols before removing missing rat orthologs:", length(unique(DEGs_df$Rat_Gene_Symbol)), "\n")

DEGs_mapped <- DEGs_df[!is.na(DEGs_df$Rat_Gene_Symbol), , drop = FALSE]

cat("Rows after requiring one mapped rat ortholog symbol:", nrow(DEGs_mapped), "\n")
cat("Unique human genes after requiring one mapped rat ortholog symbol:", length(unique(DEGs_mapped$Gene_symbol)), "\n")
cat("Unique rat ortholog symbols after requiring one mapped rat ortholog symbol:", length(unique(DEGs_mapped$Rat_Gene_Symbol)), "\n")

# True 1:1 ortholog filtering
ortholog_pairs <- unique(DEGs_mapped[, c("Gene_symbol", "Rat_Gene_Symbol")])

cat("Unique human-rat ortholog pairs:", nrow(ortholog_pairs), "\n")

human_counts <- table(ortholog_pairs$Gene_symbol)
rat_counts <- table(ortholog_pairs$Rat_Gene_Symbol)

valid_pairs <- ortholog_pairs[
  ortholog_pairs$Gene_symbol %in% names(human_counts[human_counts == 1]) &
    ortholog_pairs$Rat_Gene_Symbol %in% names(rat_counts[rat_counts == 1]),
  ,
  drop = FALSE
]

cat("Unique human genes in 1:1 ortholog pairs:", length(unique(valid_pairs$Gene_symbol)), "\n")
cat("Unique rat genes in 1:1 ortholog pairs:", length(unique(valid_pairs$Rat_Gene_Symbol)), "\n")
cat("Unique 1:1 ortholog pairs retained:", nrow(valid_pairs), "\n")

DEGs_one2one <- DEGs_mapped[
  DEGs_mapped$Gene_symbol %in% valid_pairs$Gene_symbol &
    DEGs_mapped$Rat_Gene_Symbol %in% valid_pairs$Rat_Gene_Symbol,
  ,
  drop = FALSE
]

cat("Rows after applying 1:1 ortholog filtering:", nrow(DEGs_one2one), "\n")
cat("Unique human genes after applying 1:1 ortholog filtering:", length(unique(DEGs_one2one$Gene_symbol)), "\n")
cat("Unique rat genes after applying 1:1 ortholog filtering:", length(unique(DEGs_one2one$Rat_Gene_Symbol)), "\n")

# Require rat dataset expression
DEGs_final <- DEGs_one2one[DEGs_one2one$Rat_Gene_Symbol %in% rat_genes, , drop = FALSE]

cat("Rows after requiring presence in rat dataset:", nrow(DEGs_final), "\n")
cat("Unique human genes after requiring presence in rat dataset:", length(unique(DEGs_final$Gene_symbol)), "\n")
cat("Unique rat genes after requiring presence in rat dataset:", length(unique(DEGs_final$Rat_Gene_Symbol)), "\n")

# Summary for table
initial_human_genes <- unique(DEGs_df$Gene_symbol)
mapped_human_genes <- unique(DEGs_mapped$Gene_symbol)
one2one_human_genes <- unique(DEGs_one2one$Gene_symbol)
final_human_genes <- unique(DEGs_final$Gene_symbol)

cat("\n--- Summary of human gene counts ---\n")
cat("Initial unique human genes:", length(initial_human_genes), "\n")
cat("After mapped rat ortholog filter:", length(mapped_human_genes), "\n")
cat("After true 1:1 ortholog filter:", length(one2one_human_genes), "\n")
cat("After requiring rat dataset expression:", length(final_human_genes), "\n")

cat("\n--- Genes lost at each step ---\n")
cat("Lost when requiring mapped rat ortholog:", length(setdiff(initial_human_genes, mapped_human_genes)), "\n")
cat("Lost when requiring true 1:1 ortholog:", length(setdiff(mapped_human_genes, one2one_human_genes)), "\n")
cat("Lost when requiring rat dataset expression:", length(setdiff(one2one_human_genes, final_human_genes)), "\n")

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