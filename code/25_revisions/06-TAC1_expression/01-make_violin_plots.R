rm(list = ls())
library(SpatialExperiment)
library(spatialLIBD)
library(jaffelab)
library(here)
library(sessioninfo)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(fastTopics)
library(scran)
library(scater)

dat_dir <- here("processed-data", "12_snRNA")
sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
sce <- sce[, !sce$CellType.Final == "Neuron_Ambig"]
sce <- sce[, sce$CellType.Final %in% c(
    "DRD1_MSN_A", "DRD1_MSN_B", "DRD1_MSN_C",
    "DRD1_MSN_D", "DRD2_MSN_A", "DRD2_MSN_B"
)]
load(here("processed-data","12_snRNA","070924_21colors_celltypeFinal.rda"),verbose = TRUE)

# Now make violin plots
plotDir <- here("plots", "25_revisions", "06-TAC1_expression")
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

## normalize if needed
if (!"logcounts" %in% assayNames(sce)) {
    sce <- logNormCounts(sce)
}

rownames(sce) <- rowData(sce)$gene_name
genes_to_plot <- c("TAC1", "DRD1", "DRD2")

## set factor order
sce$CellType.Final <- factor(
    sce$CellType.Final,
    levels = c(
        "DRD1_MSN_A", "DRD1_MSN_B", "DRD1_MSN_C",
        "DRD1_MSN_D", "DRD2_MSN_A", "DRD2_MSN_B"
    )
)

celltype_colors <- cluster_cols[levels(sce$CellType.Final)]

p_combined <- scater::plotExpression(
    sce,
    features = genes_to_plot,
    x = "CellType.Final",
    colour_by = "CellType.Final",
    exprs_values = "logcounts", ncol = 1
) +
    scale_color_manual(values = cluster_cols) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(x = NULL, y = "logcounts")

pdf(file.path(plotDir, "gene_expression.pdf"), width = 4, height = 6)
print(p_combined)
dev.off()