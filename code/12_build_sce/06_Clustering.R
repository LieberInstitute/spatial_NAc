# cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
# module load r_nac

#Load libraries
library(SingleCellExperiment)
library(sessioninfo)
library(HDF5Array)
library(batchelor)
library(rafalib)
library(ggplot2)
library(scater)
library(dplyr)
library(scran)
library(here)

#load sce
sce <- loadHDF5SummarizedExperiment(dir = here("processed-data",
                                               "12_snRNA",
                                               "sce_DimRed"))

print("Object loaded")
sce

#Additional sanity chekc
identical(rownames(colData(sce)),colnames(sce))
stopifnot(identical(rownames(colData(sce)),colnames(sce)))


#Comoute log counts to plot expression
print("Computing log counts")
sce <- batchelor::multiBatchNorm(sce, batch = sce$Sample)

sce

#Once more. 
identical(rownames(colData(sce)),colnames(sce))
stopifnot(identical(rownames(colData(sce)),colnames(sce)))


#Plot expression of x and y linked genes. Color by Sample
x_y_linked <- plotExpression(object = sce,
                             features = c("XIST","USP9Y","UTY"),
                             x = "Sample",
                             colour_by = "Brain_ID",
                             swap_rownames = "gene_name",
                             ncol = 1) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(plot = x_y_linked,
       filename = here("plots","12_snRNA","Expression","X_Y_genes.png"))


#Plot expression of known marker genes to be able to identify appropriate clusters
genes <- c("SNAP25","SYT1","RBFOX3", #PAN-NEURON markers
           "GAD1","GAD2","SLC32A1", #GABA markers
           "PPP1R1B","BCL11B","ISL1","FOXP2","RELN",#MSN markers
           "DRD1","PDYN","EBF1","RXFP1","TAC1","CRHR2", #D1 markers
           "DRD2","PENK","ADORA2A",#D2 markers,
           "GRM8","CHST9","OPRM1", #D1-ISLANDS
           "KIT","CHAT","SST", #INTERNEURON MARKERS
           "MOBP","MBP","OPALIN", #Oligodendrocytes
           "PDGFRA", #POLYDENDROCYTES
           "AQP4","GFAP", #Astrocytes
           "CD74") #Microglia

for(i in genes){
  print(i)
  x <- plotReducedDim(object = sce,
                      dimred = "tSNE_mnn",
                      colour_by = i,
                      swap_rownames = "gene_name") + 
  scale_color_gradientn(colours = c("lightgrey","red")) +
  ggtitle(i) +
  theme(plot.title = element_text(hjust = 0.5))
  ggsave(plot = x,
         filename = here("plots","12_snRNA","Expression","Known_Marker_Genes",
                         paste0(i,".png")))
}


snn_k_10 <- buildSNNGraph(sce, k = 10, use.dimred = "mnn",type="jaccard")

#Louvain clustering
set.seed(1234)
clust_10 <- igraph::cluster_louvain(snn_k_10,resolution=1)$membership
table(clust_10)

#Add cluster information to object
sce$k_10_louvain <- factor(clust_10)

k_10_tSNE <- plotReducedDim(object = sce,
                            dimred = "tSNE_mnn",
                            colour_by = "k_10_louvain",
                            text_by = "k_10_louvain") +
  ggtitle("k=10 louvain clustering") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_10_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_10_louvain_tSNE.png"))

#Plot tSNE by sample
for(i in unique(sce$Sample)){
  print(i)
  tSNE_sample <- plotReducedDim(object = sce[,sce$Sample == i],
                                dimred = "tSNE_mnn",
                                colour_by = "k_10_louvain") +
    ggtitle(unique(colData(sce[,sce$Sample == i])[,"Sort"])) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(plot = tSNE_sample, filename = here("plots","12_snRNA",
                                             "Dim_Red","tSNE_Sample_Specific",
                                             paste0(i,".png")))
}


#Known marker genes by cluster (violins)
for(i in genes){
  print(i)
  x <- plotExpression(object = sce,
                      features = i,
                      x = "k_10_louvain",
                      colour_by = "k_10_louvain",
                      swap_rownames = "gene_name",ncol = 1) +
    theme(legend.position = "none")
  ggsave(filename = here("plots","12_snRNA",
                         "Expression","Known_Marker_Genes","Violin_Plots",
                         paste0(i,"_violin.png")),height = 4,width = 8)
}



print("Saving sce object")
saveHDF5SummarizedExperiment(sce,
                             here("processed-data","12_snRNA",
                                  "sce_clustered"),
                             replace = TRUE)
print("Object saved")


print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
