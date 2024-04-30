#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

library(SingleCellExperiment)
library(sessioninfo)
library(HDF5Array)
library(harmony)
library(ggplot2)
library(scater)
library(scran)
library(here)

#load object
sce_neuron <- loadHDF5SummarizedExperiment(dir = here("processed-data","12_snRNA",
                                                      "sce_neuron_featureslxn"))

identical(rownames(colData(sce_neuron)),colnames(sce_neuron))

sce_neuron

#Pull the highly deviant genes
hdgs <- rownames(sce_neuron)[order(rowData(sce_neuron)$binomial_deviance, decreasing = T)][1:2000]
hdgs.symbols <- rowData(sce_neuron)$gene_name[match(hdgs, rowData(sce_neuron)$gene_id)]

print("Highly deviant genes ordered by rank")
hdgs.symbols

print("Higly deviant genes ordered alphabetically")
hdgs.symbols[order(hdgs.symbols)]

#####Dimensionality Reduction
#Remove all of the reducedDim values
reducedDim(sce_neuron,"GLMPCA_approx") <- NULL
reducedDim(sce_neuron,"tSNE") <- NULL
reducedDim(sce_neuron,"HARMONY") <- NULL
reducedDim(sce_neuron,"tSNE_HARMONY") <- NULL

#Run PCA now
print("Running PCA-subset")
Sys.time()
#Run PCA
set.seed(410)
sce_neuron <- runPCA(sce_neuron,
                     exprs_values = "binomial_pearson_residuals",
                     subset_row = hdgs,
                     ncomponents = 100,
                     name = "GLMPCA")
warnings() #Warnings seem to come from the HDF5 saved objects only. Tested this by running runPCA with a .rds object earlier
print("PCA Complete")
Sys.time()

#plotPCA by Brain
PCA_brain_neuron <- plotReducedDim(object = sce_neuron,
                                   dimred = "GLMPCA",
                                   colour_by = "Brain_ID")

ggsave(plot = PCA_brain_neuron,filename = here("plots","12_snRNA","Dim_Red","PCA_neuron_subset.png"))

#tSNE 
#t-SNE with top 50 dimensions.
print("Running tSNE")
Sys.time()
set.seed(21114)
sce_neuron <- runTSNE(sce_neuron,
                      dimred = "GLMPCA",
                      n_dimred = 50,
                      name = "tSNE")
print("tSNE complete")
Sys.time()

#plot tSNE  by Brain
tSNE_brain_neuron <- plotReducedDim(object = sce_neuron,
                                    dimred = "tSNE",
                                    colour_by = "Brain_ID")

ggsave(plot = tSNE_brain_neuron,filename = here("plots","12_snRNA","Dim_Red","tSNE_neuron_subset.png"))

#Harmony requires the PCA reduced dim to be named "PCA"
reducedDim(sce_neuron,"PCA") <- reducedDim(sce_neuron, "GLMPCA")

message("running Harmony")
sce_neuron <- RunHarmony(sce_neuron,group.by.vars = "Sample",verbose = TRUE)

#Remove the redundant PCA reducedDim
reducedDim(sce_neuron, "PCA") <- NULL

#tSNE post-mnn with top 50 dimensions
print("Running tSNE post-Harmony")
set.seed(21114)
Sys.time()
sce_neuron <- runTSNE(sce_neuron,
                      dimred = "HARMONY",
                      n_dimred = 50,
                      name = "tSNE_HARMONY")
Sys.time()

#plot tSNE  by Brain
tSNE_brain_neuron_harmony <- plotReducedDim(object = sce_neuron,
                                            dimred = "tSNE_HARMONY",
                                            colour_by = "Brain_ID")

ggsave(plot = tSNE_brain_neuron_harmony,filename = here("plots","12_snRNA","Dim_Red","tSNE_neuron_subset_harmony.png"))

#Save object
print("Saving object")
Sys.time()
saveHDF5SummarizedExperiment(sce_neuron,
                             here("processed-data",
                                  "12_snRNA",
                                  "sce_neuron_DimRed"),
                             replace = TRUE)

#########
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
