library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
library(scater)
library(scran)
library(SingleCellExperiment)
library(RcppML)
library(viridis)
library(pals)
set.seed(1234)

# Read data and create Seurat object
dat_dir <- here::here("processed-data", "12_snRNA")
gene_selection_strategy <- "top5000_HVGs"
res_dir <- here::here("processed-data", "19_snRNAseq_NMF", "RcppML", gene_selection_strategy)
plot_dir <- here::here("plots", "19_snRNAseq_NMF", "RcppML", gene_selection_strategy)
dir.create(res_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
# Convert rownames from ensembl ID to symbol
sce$CellType.Final[sce$CellType.Final == "T-Cell"] <- "T_cell"
rownames(sce) <- rowData(sce)$gene_name
rownames(sce) <- make.names(rownames(sce), unique = TRUE)
# Get counts and metadata that are required to initialize Seurat object
counts <- counts(sce)
metadata <- colData(sce)

# Initialize the Seurat object
sobj <- CreateSeuratObject(counts = counts, project = "NAc_snRNAseq", min.cells = 3, min.features = 200)
metadata <- metadata[rownames(metadata) %in% colnames(sobj), ]
metadata <- metadata[match(colnames(sobj), rownames(metadata)), ]
sobj$Sample <- metadata$Sample
sobj$Brain_ID <- metadata$Brain_ID
sobj$percent_mito <- metadata$subsets_Mito_percent
sobj$CellType <- metadata$CellType.Final
sobj$Sort <- metadata$Sort
sobj$snRNA_date <- metadata$snRNA_data
sobj$Chromium_cDNA_date <- metadata$Chromium_cDNA_date
sobj$Chromium_library_date <- metadata$Chromium_library_date
sobj$Chromium_kit <- metadata$Chromium_kit
sobj$Library_Kit <- metadata$Library_Kit

# Process single cell data
sobj <- sobj %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)

sobj <- sobj %>%
  FindNeighbors(dims = 1:20, verbose = FALSE) %>%
  FindClusters(resolution = 1.0, verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)

if(gene_selection_strategy == "all_genes"){
  logNorm_mat <- as(sobj[["RNA"]]$data ,"sparseMatrix")
}else{
  logNorm_mat <- as(sobj[["RNA"]]$data[VariableFeatures(sobj), ] ,"sparseMatrix")
}

K <- 30
model <- nmf(logNorm_mat, k = K, seed = 1:10, L1 = c(0.1, 0.1))
saveRDS(model, file.path(res_dir, "nmf_model.rds"))

summary_df <- summary(model, group_by = sobj$Brain_ID, stat = "mean")
summary_df$factor <- factor(summary_df$factor, levels = paste0("nmf", c(1:K)))
colnames(summary_df)[1] <- "Sample_ID"
pdf(file.path(plot_dir, "weights_by_sampleID.pdf"), width = 12, height = 5)
ggplot(summary_df, aes(x = factor, y = stat, fill = Sample_ID)) + geom_bar(stat = "identity", position = "fill") + 
scale_fill_manual(values = pal_jco()(10))  + xlab("") +
theme_pubr() + labs_pubr() + ylab("Mean NMF weight") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()

summary_df <- summary(model, group_by = sobj$CellType, stat = "mean")
summary_df$factor <- factor(summary_df$factor, levels = paste0("nmf", c(1:K)))
colnames(summary_df)[1] <- "Cell_type"
pdf(file.path(plot_dir, "weights_by_cellType.pdf"), width = 12, height = 5)
ggplot(summary_df, aes(x = factor, y = stat, fill = Cell_type)) + geom_bar(stat = "identity", position = "fill") + 
scale_fill_manual(values = cols25(n = 25))  + xlab("") +
theme_pubr() + labs_pubr() + ylab("Mean NMF weight") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()

# Add an assay to the single cell data
stopifnot(all.equal(colnames(sobj), colnames(model$h)))
sobj[["RcppML"]] <- CreateAssayObject(counts = model$h)
color_palette <- viridis(n=10)
pattern_names <- rownames(model$h)

pdf(file.path(plot_dir, "UMAP_sampleID.pdf"), width = 6, height = 5)
DimPlot(sobj, reduction = "umap", group.by = "Sample")
dev.off()

pdf(file.path(plot_dir, "UMAP_BrainID.pdf"), width = 6, height = 5)
DimPlot(sobj, reduction = "umap", group.by = "Brain_ID", cols = pal_jco()(10))
dev.off()

pdf(file.path(plot_dir, "UMAP_cellType.pdf"), width = 10, height = 5)
DimPlot(sobj, reduction = "umap", group.by = "CellType", cols = cols25(n = 25))
dev.off()

plot_umap_list <- list()
pindex <- 1
nPlot_size <- 6
start_index <- 1
stop_index <- start_index + nPlot_size - 1
while(start_index <= K){
    plot_umap_list[[pindex]] <- FeaturePlot(sobj, pattern_names[start_index:stop_index], cols=color_palette, reduction = "umap", ncol = 3) & NoLegend()
    pindex <- pindex + 1
    start_index <- stop_index + 1
    stop_index <- start_index + nPlot_size - 1
    stop_index <- min(stop_index, K)
}
pdf(file.path(plot_dir, "UMAP_NMF_factors.pdf"), width = 12, height = 8)
for(i in c(1:length(plot_umap_list))){
    print(plot_umap_list[[i]])
}
dev.off()
