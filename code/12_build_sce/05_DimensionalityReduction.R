#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(harmony)
library(scater)
library(scran)
library(scry)
library(here)


#Load sce 
sce <- readRDS(here("processed-data","12_snRNA","sce_featureselection.Rds"))

sce

#Take top 4000 highly deviant genes
hdgs <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing = T)][1:4000]
hdgs.symbols <- rowData(sce)$gene_name[match(hdgs, rowData(sce)$gene_id)]

#Check if Dopamine receptors and PPP1R1B are among the highly deviant genes
c("DRD1","DRD2","PPP1R1B") %in% hdgs.symbols

print("Running PCA")
#Run PCA
set.seed(1010)
sce <- runPCA(sce,
              exprs_values = "binomial_pearson_residuals",
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
print("PCA Complete")

#t-SNE with 100 dimensions.
print("Running tSNE") 
set.seed(1010)
sce <- runTSNE(sce,
               dimred = "GLMPCA_approx",
               n_dimred = 100,
               name = "tSNE")
print("tSNE complete")

#Visualize t-SNE with points colored by sample
tsne_sample <- plotReducedDim(sce,
                              dimred = "tSNE", 
                              colour_by = "Sample",
                              point_alpha = 0.3)
ggsave(tsne_sample,filename = here("plots","12_snRNA","Dim_Red","tSNE_Sample_noCorrection.png"))
#Looks like some sections af the tSNE are dominated by just a few of the samples. 

#Visualize t-SNE with points colored by sample
tsne_brainid <- plotReducedDim(sce,
                              dimred = "tSNE",
                              colour_by = "Brain_ID",
                              point_alpha = 0.3)
ggsave(tsne_brainid,filename = here("plots","12_snRNA","Dim_Red","tSNE_BrainID_noCorrection.png"))

#Plot by snRNA_date
#Misspelled as snRNA_data, will fix later. 
tsne_snRNA_date <- plotReducedDim(sce,
                                  dimred = "tSNE", 
                                  colour_by = "snRNA_data",
                                  point_alpha = 0.3)
ggsave(tsne_snRNA_date,filename = here("plots","12_snRNA","Dim_Red","tSNE_snRNA_date_noCorrection.png"))

#Plot by sort type
tsne_sort <- plotReducedDim(sce,
                            dimred = "tSNE", 
                            colour_by = "Sort",
                            point_alpha = 0.3)
ggsave(tsne_sort,filename = here("plots","12_snRNA","Dim_Red","tSNE_SortType_noCorrection.png"))

#UMAP  with 100 dimensions
print("Running UMAP")
set.seed(1010)
sce <- runUMAP(sce,
               dimred = "GLMPCA_approx",
               n_dimred = 100,
               name = "umap")

#UMAP by Sample 
umap_Sample <- plotReducedDim(sce,
                              dimred = "umap",
                              colour_by = "Sample",
                              point_alpha = 0.3)
ggsave(umap_Sample,
       filename = here("plots","12_snRNA","Dim_Red","umap_Sample_noCorrection.png"))

#UMAP by Brain_ID 
umap_BrainID <- plotReducedDim(sce,
                               dimred = "umap",
                               colour_by = "Brain_ID",
                               point_alpha = 0.3)
ggsave(umap_BrainID,
       filename = here("plots","12_snRNA","Dim_Red","umap_BrainID_noCorrection.png"))

#UMAP by snRNA_date
umap_snRNA_date <- plotReducedDim(sce,
                                  dimred = "umap",
                                  colour_by = "snRNA_data",
                                  point_alpha = 0.3)
ggsave(umap_snRNA_date,filename = here("plots","12_snRNA","Dim_Red","umap_snRNA_date_noCorrection.png"))

#UMAP by sort type
umap_sort <- plotReducedDim(sce,
                            dimred = "umap",
                            colour_by = "Sort",
                            point_alpha = 0.3)
ggsave(umap_sort,filename = here("plots","12_snRNA","Dim_Red","umap_SortType_noCorrection.png"))

print("umap complete")

#########HARMONY
#Harmony requires the PCA reduced dim to be named "PCA"
reducedDim(sce,"PCA") <- reducedDim(sce, "GLMPCA_approx")

message("running Harmony")
set.seed(1010)
sce <- RunHarmony(sce,group.by.vars = "Sample",verbose = TRUE)

#Remove the redundant PCA reducedDim
reducedDim(sce, "PCA") <- NULL

print("Dimensions of HARMONY matrix")
dim(reducedDim(sce, "HARMONY"))

#tSNE post-HARMONY with top 100 dimensions (or all that were calculated during PCA step)
print("Running tSNE post-Harmony")
set.seed(1010)
sce <- runTSNE(sce,
               dimred = "HARMONY",
               n_dimred = 100,
               name = "tSNE_HARMONY")

#Corrected tSNE by Sample 
tSNE_HARMONY_Sample <- plotReducedDim(sce,
                                      dimred = "tSNE_HARMONY", 
                                      colour_by = "Sample",
                                      point_alpha = 0.3)
ggsave(tSNE_HARMONY_Sample,
       filename = here("plots","12_snRNA","Dim_Red","tSNE_Sample_HARMONY.png"))

#Corrected tSNE by Brain_ID 
tSNE_HARMONY_BrainID <- plotReducedDim(sce,
                                       dimred = "tSNE_HARMONY",
                                       colour_by = "Brain_ID",
                                       point_alpha = 0.3)
ggsave(tSNE_HARMONY_BrainID,
       filename = here("plots","12_snRNA","Dim_Red","tSNE_BrainID_HARMONY.png"))

#Plot by snRNA_date
#Misspelled as snRNA_data, will fix later. 
tsne_HARMONY_snRNA_date <- plotReducedDim(sce,
                                          dimred = "tSNE_HARMONY",
                                          colour_by = "snRNA_data",
                                          point_alpha = 0.3)
ggsave(tsne_HARMONY_snRNA_date,filename = here("plots","12_snRNA","Dim_Red","tSNE_snRNA_date_HARMONY.png"))

#Plot by sort type
tsne_HARMONY_sort <- plotReducedDim(sce,
                            dimred = "tSNE_HARMONY",
                            colour_by = "Sort",
                            point_alpha = 0.3)
ggsave(tsne_HARMONY_sort,filename = here("plots","12_snRNA","Dim_Red","tSNE_SortType_HARMONY.png"))

print("tSNE_HARMONY complete")

#UMAP post Harmony with 100 dimensions
print("Running UMAP post-HARMONY")
set.seed(1010)
sce <- runUMAP(sce,
               dimred = "HARMONY",
               n_dimred = 100,
               name = "umap_HARMONY")

#Corrected tSNE by Sample 
umap_HARMONY_Sample <- plotReducedDim(sce,
                                      dimred = "umap_HARMONY",
                                      colour_by = "Sample",
                                      point_alpha = 0.3)
ggsave(umap_HARMONY_Sample,
       filename = here("plots","12_snRNA","Dim_Red","umap_Sample_HARMONY.png"))

#Corrected tSNE by Brain_ID 
umap_HARMONY_BrainID <- plotReducedDim(sce,
                                       dimred = "umap_HARMONY",
                                       colour_by = "Brain_ID",
                                       point_alpha = 0.3)
ggsave(umap_HARMONY_BrainID,
       filename = here("plots","12_snRNA","Dim_Red","umap_BrainID_HARMONY.png"))

#Plot by snRNA_date
#Misspelled as snRNA_data, will fix later. 
umap_HARMONY_snRNA_date <- plotReducedDim(sce,
                                          dimred = "umap_HARMONY",
                                          colour_by = "snRNA_data",
                                          point_alpha = 0.3)
ggsave(umap_HARMONY_snRNA_date,filename = here("plots","12_snRNA","Dim_Red","umap_snRNA_date_HARMONY.png"))

#Plot by sort type
umap_HARMONY_sort <- plotReducedDim(sce,
                            dimred = "umap_HARMONY",
                            colour_by = "Sort",
                            point_alpha = 0.3)
ggsave(umap_HARMONY_sort,filename = here("plots","12_snRNA","Dim_Red","umap_SortType_HARMONY.png"))

print("umap_HARMONY complete")

#Save object with Dimensionality Reduction
#Switching to save with HDF5 here because at this point there is so much information within the sce
print("Saving sce object")
saveRDS(sce,
        file = here("processed-data","12_snRNA","sce_DimRed.Rds"))

###Reproduciblity
print("Object saved")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
