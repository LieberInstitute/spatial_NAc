library(SpatialExperiment)
library(sessioninfo)
library(ggplot2)
library(Seurat)
library(escheR)
library(here)

#Load spe object containing normalized coutns
spe <- readRDS(here("processed-data","xenium","SPEs","spe_NormCounts.Rds"))

#set logcounts
logcounts(spe) <- assay(spe,"nucleus_normcounts")

#Convert to seurat
seurat_spe <- as.Seurat(spe,counts = "counts",data = "logcounts")

Assays(seurat_spe)

#Run PCA 
seurat_spe <- FindVariableFeatures(seurat_spe)
seurat_spe <- ScaleData(seurat_spe)
seurat_spe <- RunPCA(seurat_spe,seed.use = 1314,npcs = 50,reduction.name = "seurat_pca")
seurat_spe <- RunUMAP(seurat_spe,seed.use = 1314, dims= 1:50,reduction = "seurat_pca") 

#Save the object
saveRDS(seurat_spe, here("processed-data", "xenium", 
                          "SPEs", "seurat_xenium.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
