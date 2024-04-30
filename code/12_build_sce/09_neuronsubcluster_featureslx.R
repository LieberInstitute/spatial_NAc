#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac


library(SingleCellExperiment)
library(sessioninfo)
library(HDF5Array)
library(bluster)
library(ggplot2)
library(scater)
library(scran)
library(here)
library(scry)

##############
#Load sce
sce <- loadHDF5SummarizedExperiment(dir = here("processed-data","12_snRNA",
                                               "sce_clustered"))

dim(sce)


identical(rownames(colData(sce)),colnames(sce))


sce
#########################

###Annotation
annotation_df <- data.frame(cluster= 1:17,
                            celltype = c("Oligo","DRD1_MSN_A","Microglia","Polydendrocyte",
                                         "DRD2_MSN_A","Ependymal","Astrocyte","DRD1_MSN_B",
                                         "Endothelial","PVALB_Inh","DRD2_MSN_B","DRD1_MSN_C",
                                         "VIP_CCK_Inh","GLP1R_Inh","GABA_Glut","SST_Inh",
                                         "T_Cells"))

#add celltype info
sce$CellType.Broad <- annotation_df$celltype[match(sce$k_10_louvain_pt25,
                                                   annotation_df$cluster)]
#factorize. 
sce$CellType.Broad <- factor(sce$CellType.Broad,
                             levels = c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C",
                                        "DRD2_MSN_A","DRD2_MSN_B",
                                        "GLP1R_Inh","PVALB_Inh","VIP_CCK_Inh","SST_Inh","GABA_Glut",
                                        "Astrocyte","Ependymal","Oligo","Polydendrocyte","Microglia",
                                        "Endothelial","T_Cells"))

#Set up new cluster colors. 
cluster_cols <- Polychrome::createPalette(length(unique(sce$CellType.Broad)),
                                          c("#D81B60", "#1E88E5","#FFC107","#009E73"))
names(cluster_cols) <- unique(sce$CellType.Broad)

save(cluster_cols,file = here("processed-data","12_snRNA","CellType_Broad_cluster_cols.rda"))

#Annotate the tSNE
tsne_broad <- plotReducedDim(object = sce,
                             dimred = "tSNE_HARMONY",
                             colour_by = "CellType.Broad",
                             text_by = "CellType.Broad") +
    scale_color_manual(values = cluster_cols)

ggsave(filename = here("plots","12_snRNA","Dim_Red","tSNE_CellType_Broad.png"),plot = tsne_broad)

#Make a colData column for neurons or glia
sce$Neuron_NonNeuron <- ifelse(sce$CellType.Broad %in% c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C",
                                                         "DRD2_MSN_A","DRD2_MSN_B","GLP1R_Inh",
                                                         "PVALB_Inh","VIP_CCK_Inh","SST_Inh",
                                                         "GABA_Glut"),
                               "Neuron",
                               "NonNeuron")

table(sce$Neuron_NonNeuron)


###
#Subset neurons
sce_neuron <- sce[, sce$Neuron_NonNeuron == "Neuron"]

sce_neuron

rm(sce)

####Feature selection
# feature selection
print("Running devianceFeatureSelection")
Sys.time()
set.seed(410)
sce_neuron <- devianceFeatureSelection(sce_neuron,
                                       assay = "counts",
                                       fam = "binomial",
                                       sorted = FALSE,
                                       batch = sce_neuron$Sample)
print("deviance FeatureSelection complete")

#Plot binomial deviance
pdf(here("plots","12_snRNA","featureSelxn_binomialDeviance-byGene_neuronsubset.pdf"))
plot(sort(rowData(sce_neuron)$binomial_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylab = "binomial deviance"
)
abline(v = 2000,lty = 2, col = "red")
dev.off()

#Calculate null residuals
#Pearson residuals. 
print("Running nullResiduals")
Sys.time()
sce_neuron <- nullResiduals(sce_neuron,
                            assay = "counts", 
                            fam   = "binomial", 
                            type  = "pearson")
print("nullResiduals complete")
Sys.time()

###Save the object and the clustersweep output
print("Saving Objects")
Sys.time()
saveHDF5SummarizedExperiment(sce_neuron,
                             here("processed-data",
				  "12_snRNA",
                                  "sce_neuron_featureslxn"),
                             replace = TRUE)
warnings()
#########
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
