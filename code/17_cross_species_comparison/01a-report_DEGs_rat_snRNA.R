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

#spec <- matrix(
#  c(
#    "subset_neurons", "n", 1, "logical", "Perform correlation only for neuronal cell types?"
#  ),
#  byrow = TRUE, ncol = 5
#)
#opt <- getopt(spec)
opt <- list()
opt$subset_neurons <- TRUE

if (is.null(opt$subset_neurons)) {
  opt$subset_neurons <- FALSE
}

dat_dir <- here("processed-data", "12_snRNA")
reg_dir <- here("processed-data", "10_post_clustering_analysis", "02_spatial_registration_sn")

# -----------------------------
# Read data
# -----------------------------
sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
sce <- sce[, !sce$CellType.Final == "Neuron_Ambig"]

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

sn_registration <- readRDS(file.path(reg_dir, "sn_cellType_registration.rds"))
sn_enrichment <- sn_registration[["enrichment"]]

if (opt$subset_neurons) {
  select_cell_types <- c(
    "DRD1_MSN_A", "DRD1_MSN_B", "DRD1_MSN_C", "DRD1_MSN_D", "DRD2_MSN_A",
    "DRD2_MSN_B", "Excitatory", "Inh_A", "Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F"
  )
  rat_genes <- rownames(sce_ext_neurons)
} else {
  select_cell_types <- c(
    "DRD1_MSN_A", "DRD1_MSN_B", "DRD1_MSN_C", "DRD1_MSN_D", "DRD2_MSN_A",
    "DRD2_MSN_B", "Excitatory", "Inh_A", "Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F",
    "Astrocyte_A", "Astrocyte_B", "Microglia", "OPC", "Oligo", "Ependymal", "Endothelial"
  )
  rat_genes <- rownames(sce_ext)
}

sce <- sce[, sce$CellType.Final %in% select_cell_types]

# -----------------------------
# Build DEG table
# -----------------------------
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

# -----------------------------
# Add human gene symbols
# -----------------------------
geneData <- rowData(sce)
DEGs_df$Gene_symbol <- NA_character_

for (i in seq_len(nrow(DEGs_df))) {
  DEGs_df$Gene_symbol[i] <- geneData$gene_name[geneData$gene_id == DEGs_df$Gene_ID[i]][1]
}

cat("Unique human genes before removing missing symbols:", length(unique(DEGs_df$Gene_symbol)), "\n")

DEGs_df <- DEGs_df[!is.na(DEGs_df$Gene_symbol), , drop = FALSE]

cat("Rows after removing missing human symbols:", nrow(DEGs_df), "\n")
cat("Unique human genes after removing missing symbols:", length(unique(DEGs_df$Gene_symbol)), "\n")

# -----------------------------
# Add rat ortholog symbols
# -----------------------------
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

# -----------------------------
# Unique ortholog pairs
# -----------------------------
ortholog_pairs <- unique(DEGs_mapped[, c("Gene_symbol", "Rat_Gene_Symbol")])

cat("Unique human-rat ortholog pairs:", nrow(ortholog_pairs), "\n")

# -----------------------------
# True 1:1 ortholog filtering
# -----------------------------
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

# -----------------------------
# Require rat expression
# -----------------------------
DEGs_final <- DEGs_one2one[DEGs_one2one$Rat_Gene_Symbol %in% rat_genes, , drop = FALSE]

cat("Rows after requiring presence in rat dataset:", nrow(DEGs_final), "\n")
cat("Unique human genes after requiring presence in rat dataset:", length(unique(DEGs_final$Gene_symbol)), "\n")
cat("Unique rat genes after requiring presence in rat dataset:", length(unique(DEGs_final$Rat_Gene_Symbol)), "\n")

# -----------------------------
# Optional summary of losses
# -----------------------------
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