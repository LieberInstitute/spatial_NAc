####snRNA-seq NMF pattern projection to Visium
library(tidyverse)
library(RcppML)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(here)
library(scater)
library(scran)
library(BiocParallel)
library(BiocSingular)
library(spatialNAcUtils)
library(jaffelab)
library(projectR)
library(scater)
library(scran)
library(dittoSeq)
library(escheR)
library(getopt)
library(Seurat)

opt <- list()
opt$data <- "rat_case_control_morphine_repeated"
print(opt$data)
res_dir <- here("processed-data", "16_transfer_learning", "02_target_projections", "projection_of_human_factors")
plot_dir <- here("plots", "16_transfer_learning", "02_target_projections", "projection_human_factors", opt$data)

sce <- readRDS(file = file.path(res_dir, paste0(opt$data, ".rds")))

pdf(file.path(plot_dir, "violin_plots.pdf"), width = 9, height = 3)
VlnPlot(sce, group.by = "CellType", features = "nmf44", pt.size = 0) + theme(legend.position = "none")
VlnPlot(sce, group.by = "CellType", features = "nmf34", pt.size = 0) + theme(legend.position = "none")
VlnPlot(sce, group.by = "CellType", features = "nmf35", pt.size = 0) + theme(legend.position = "none")
VlnPlot(sce, group.by = "CellType", features = "nmf3", pt.size = 0) + theme(legend.position = "none")
VlnPlot(sce, group.by = "CellType", features = "nmf4", pt.size = 0) + theme(legend.position = "none")
VlnPlot(sce, group.by = "CellType", features = "nmf7", pt.size = 0) + theme(legend.position = "none")
VlnPlot(sce, group.by = "CellType", features = "nmf10", pt.size = 0) + theme(legend.position = "none")
dev.off()

# Examine cell types that are selected based on the NMF mappings
nmf10_summ <- table(sce$CellType[sce$nmf10 > 3e-5], sce$treatment[sce$nmf10 > 3e-5])
nmf10_summ <- reshape2::melt(nmf10_summ)
pdf(file.path(plot_dir, "selected_cells_nmf10.pdf"), width = 4, height = 7)
ggplot(nmf10_summ, aes(x = Var1, y = value, fill = Var2)) + geom_bar(position = "dodge", stat = "identity") + ggtitle("nmf10") + 
coord_flip() + theme_minimal() + xlab("Cell Type") + ylab("Frequency") + guides(fill=guide_legend(title="Treatment"))
dev.off()

nmf7_summ <- table(sce$CellType[sce$nmf7 > 3e-5], sce$treatment[sce$nmf7 > 3e-5])
nmf7_summ <- reshape2::melt(nmf7_summ)
pdf(file.path(plot_dir, "selected_cells_nmf7.pdf"), width = 4, height = 7)
ggplot(nmf7_summ, aes(x = Var1, y = value, fill = Var2)) + geom_bar(position = "dodge", stat = "identity") + ggtitle("nmf7") + 
coord_flip() + theme_minimal() + xlab("Cell Type") + ylab("Frequency") + guides(fill=guide_legend(title="Treatment"))
dev.off()

nmf44_summ <- table(sce$CellType[sce$nmf44 > 5e-5], sce$treatment[sce$nmf44 > 5e-5])
nmf44_summ <- reshape2::melt(nmf44_summ)
pdf(file.path(plot_dir, "selected_cells_nmf44.pdf"), width = 4, height = 7)
ggplot(nmf44_summ, aes(x = Var1, y = value, fill = Var2)) + geom_bar(position = "dodge", stat = "identity") + ggtitle("nmf44") + 
coord_flip() + theme_minimal() + xlab("Cell Type") + ylab("Frequency") + guides(fill=guide_legend(title="Treatment"))
dev.off()

nmf_balance_DE <- function(sce, nmf_col = "nmf10", threshold = 3e-5, 
                           conds = c("morphine","saline"),
                           plot_dir = ".",
                           seed = 123) {
  # Dependencies
  stopifnot(inherits(sce, "Seurat"))   # assumes you’re using Seurat object
  require(Seurat)
  require(tidyverse)
  
  set.seed(seed)
  
  # 1) Select cells above threshold
  meta <- sce@meta.data
  meta$cell <- rownames(meta)
  meta <- meta[meta$treatment %in% conds & meta[[nmf_col]] > threshold, , drop = FALSE]
  
  balance_one <- function(df_ct) {
    m_idx <- which(df_ct$treatment == conds[1])
    s_idx <- which(df_ct$treatment == conds[2])
    k <- min(length(m_idx), length(s_idx))
    if (k == 0) return(df_ct[integer(0), , drop = FALSE])
    m_take <- if (length(m_idx) > k) sample(m_idx, k) else m_idx
    s_take <- if (length(s_idx) > k) sample(s_idx, k) else s_idx
    df_ct[c(m_take, s_take), , drop = FALSE]
  }
  # split by CellType, balance, and row-bind
  split_ct <- split(meta, meta$CellType, drop = TRUE)
  balanced_meta_list <- lapply(split_ct, balance_one)
  balanced_meta <- do.call(rbind, balanced_meta_list)

  if (nrow(balanced_meta) == 0) {
    stop("After balancing, zero cells remain (no CellType present in both conditions above threshold).")
  }
  balanced_cells <- balanced_meta$cell

  # (Optional) text summary instead of printing a tibble (avoids your error)
  tab <- with(balanced_meta, table(CellType, treatment))
  print(tab)
  
  # 3) DE test
  Idents(sce) <- sce$treatment
  markers <- FindMarkers(
    object = sce,
    ident.1 = conds[1],
    ident.2 = conds[2],
    subset.cells = balanced_cells,
    group.by = "treatment",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.1
  )
  
  markers$gene <- rownames(markers)
  markers$neglog10_padj <- -log10(markers$p_val_adj + 1e-300)
  markers$sig <- with(markers, p_val_adj < 0.05 & abs(avg_log2FC) >= 0.25)

  saveDir <- here("processed-data", "16_transfer_learning", "05_human_projection_DEGs", opt$data)
  plotDir2 <- here("plots", "16_transfer_learning", "05_human_projection_DEGs", opt$data)
  dir.create(saveDir)  
  dir.create(plotDir2, recursive = TRUE)
  # Save DE results
  de_path <- file.path(saveDir, paste0("DE_balanced_", nmf_col, "_", conds[1], "_vs_", conds[2], ".tsv"))
  write.table(markers, de_path, sep="\t", quote=FALSE, row.names=FALSE)
  message("DE table saved to: ", de_path)
  
  # Volcano plot
  top_genes <- markers %>% arrange(p_val_adj) %>% slice_head(n=15)
  volcano <- EnhancedVolcano(markers, lab = markers$gene, x = 'avg_log2FC', y = 'p_val_adj', FCcutoff = 0.5)
  
  pdf_path <- file.path(plotDir2, paste0("volcano_balanced_", nmf_col, "_", conds[1], "_vs_", conds[2], ".pdf"))
  pdf(pdf_path, width=6, height=5)
  print(volcano)
  dev.off()
  message("Volcano plot saved to: ", pdf_path)
  
  return(markers)
}
