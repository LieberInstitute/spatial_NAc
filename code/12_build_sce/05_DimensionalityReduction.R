#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

library(SingleCellExperiment)
library(sessioninfo)
library(HDF5Array)
library(ggplot2)
library(scater)
library(scran)
library(scry)
library(here)

#Load sce 
sce <- loadHDF5SummarizedExperiment(dir = here("processed-data","12_snRNA",
                                               "sce_featureselection"))

sce

#Take top 2000 highly deviant genes
hdgs <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing = T)][1:2000]
hdgs.symbols <- rowData(sce)$gene_name[match(hdgs, rowData(sce)$gene_id)]

#Checkl if Dopamine receptors and PPP1R1B are among the highly deviant genes
c("DRD1","DRD2","PPP1R1B") %in% hdgs.symbols

#Run PCA
sce <- runPCA(sce,
              exprs_values = "binomial_deviance_residuals",
              subset_row = hdgs, 
              ncomponents = 100,
              name = "GLMPCA_approx")

#PCA plot of top 6 PCs
PCA_plots <- plotReducedDim(sce,
                            dimred = "GLMPCA_approx", 
                            colour_by = "Sample",
                            ncomponents = 6, 
                            point_alpha = 0.3)
ggsave(PCA_plots,filename = here("plots","12_snRNA","Dim_Red","Top6_PCAs.png"))


#t-SNE with top 50 dimensions. 
set.seed(100111001)
sce <- runTSNE(sce,
               dimred = "GLMPCA_approx",
               n_dimred = 50,
               name = "tSNE")

#Visualize t-SNE with points colored by sample
tsne_sample <- plotReducedDim(sce,
                              dimred = "tSNE", 
                              colour_by = "Sample",
                              point_alpha = 0.3)
ggsave(tsne_sample,filename = here("plots","12_snRNA","Dim_Red","tSNE_Sample_noMNN.png"))
#Looks like some sections af the tSNE are dominated by just a few of the samples. 

#Plot by snRNA_date
tsne_snRNA_date <- plotReducedDim(sce,
                                  dimred = "tSNE", 
                                  colour_by = "snRNA_date",
                                  point_alpha = 0.3)
ggsave(tsne_snRNA_date,filename = here("plots","12_snRNA","Dim_Red","tSNE_snRNA_date_noMNN.png"))

#Plot by sort type
tsne_sort <- plotReducedDim(sce,
                            dimred = "tSNE", 
                            colour_by = "Sort",
                            point_alpha = 0.3)
ggsave(tsne_sort,filename = here("plots","12_snRNA","Dim_Red","tSNE_SortType_noMNN.png"))

#Batch correction most likely needed here. 
#Run batch correction with mutual nearest neighbors. 
glmpca_mnn <- batchelor::reducedMNN(reducedDim(sce, "GLMPCA_approx"),
                                    batch=as.factor(sce$Sample))

#Add mnn to the object
reducedDim(sce,"mnn") <- glmpca_mnn$corrected

#tSNE post-mnn with top 50 dimensions
set.seed(100111001)
sce <- runTSNE(sce,
               dimred = "mnn",
               n_dimred = 50,
               name = "tSNE_mnn")

tSNE_mnn_Sample <- plotReducedDim(sce,
                                  dimred = "tSNE_mnn", 
                                  colour_by = "Sample",
                                  point_alpha = 0.3)
ggsave(tSNE_mnn_Sample,
       filename = here("plots","12_snRNA","Dim_Red","tSNE_Sample_mnn.png"))

#Save object with Dimensionality Reduction
saveHDF5SummarizedExperiment(sce,
                             here("processed-data","12_snRNA",
                                  "sce_DimRed"),
                             replace = TRUE)

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
