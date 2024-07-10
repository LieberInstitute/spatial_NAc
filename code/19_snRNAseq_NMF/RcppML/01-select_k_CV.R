library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(scater)
library(scran)
library(SingleCellExperiment)
library(RcppML)
set.seed(1234)

# Read data and create Seurat object
dat_dir <- here::here("processed-data", "12_snRNA")
res_dir <- here::here("processed-data", "19_snRNAseq_NMF", "RcppML")
plot_dir <- here::here("plots", "19_snRNAseq_NMF", "RcppML")

sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
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
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)

pdf(file.path(plot_dir, "snRNA_seq_PCs.pdf"), height = 6, width = 10)
# Check Elbow plot to determine relevant number of PCs
ElbowPlot(sobj, ndims = 50)
DimPlot(sobj, reduction = "pca", group.by = "Sample")
DimPlot(sobj, reduction = "pca", group.by = "Brain_ID")
DimPlot(sobj, reduction = "pca", group.by = "Sort")
DimPlot(sobj, reduction = "pca", group.by = "snRNA_date")
FeaturePlot(sobj, reduction = "pca", features = c("nCount_RNA"))
FeaturePlot(sobj, reduction = "pca", features = c("nFeature_RNA"))
FeaturePlot(sobj, reduction = "pca", features = c("percent_mito"))
dev.off()

# Get the proportion of variance of each PC explained by meta-data
pca_embeddings <- Embeddings(sobj, reduction = "pca")[, 1:20]
pca_res <- matrix(NA, nrow = dim(pca_embeddings)[2], ncol = 11)
for(i in c(1:dim(pca_embeddings)[2])){
  cat(i, "\n")
  pca_res[i, 1] <- summary(lm(pca_embeddings[ ,i] ~ sobj$nCount_RNA))$adj.r.squared
  pca_res[i, 2] <- summary(lm(pca_embeddings[ ,i] ~ sobj$nFeature_RNA))$adj.r.squared
  pca_res[i, 3] <- summary(lm(pca_embeddings[ ,i] ~ sobj$percent_mito))$adj.r.squared
  pca_res[i, 4] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$Sample)))$adj.r.squared
  pca_res[i, 5] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$Brain_ID)))$adj.r.squared
  pca_res[i, 6] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$CellType)))$adj.r.squared
  pca_res[i, 7] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$Sort)))$adj.r.squared
  pca_res[i, 8] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$snRNA_date)))$adj.r.squared
  pca_res[i, 9] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$Chromium_cDNA_date)))$adj.r.squared
  pca_res[i, 10] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$Chromium_library_date)))$adj.r.squared
  pca_res[i, 11] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$Chromium_kit)))$adj.r.squared
}
colnames(pca_res) <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "Sample", "Brain_ID", "Cell_Type", "Sort", "snRNA_date", "Chromium_cDNA_date", "Chromium_library_date", "Chromium_kit")
rownames(pca_res) <- paste0("PC_", c(1:20))
pca_res <- reshape2::melt(pca_res)

pdf(file.path(plot_dir, "snRNA_seq_PC_drivers.pdf"), height = 6, width = 10)
ggplot(pca_res, aes(x = Var1, y = Var2, fill = value)) + geom_tile(color = "black") + 
scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) + coord_fixed() + xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()

sobj <- sobj %>%
  FindNeighbors(dims = 1:20, verbose = FALSE) %>%
  FindClusters(resolution = 1.0, verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)

# Determine the optimal rank for NMF
logNorm_mat <- as(sobj[["RNA"]]$data[VariableFeatures(sobj), ] ,"sparseMatrix")
cv <- crossValidate(logNorm_mat, k = 3:50, method = "robust", reps = 8, seed = 123)

pdf(file.path(plot_dir, "robust_cross_validation.pdf"), width = 8, height = 5)
plot(cv) + ggtitle("bi-cross-validation (Robust)")
dev.off()

saveRDS(cv, file.path(res_dir, "robust_cv.rds"))
cv <- data.frame(cv)
ggplot(cv, aes(x = k, y = value, color = rep)) + geom_point() + geom_line() + theme_pubr() + facet_wrap(~rep) + ylim(c(0, 2))