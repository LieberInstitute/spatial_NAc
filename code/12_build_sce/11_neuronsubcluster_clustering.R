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
sce_neuron <- loadHDF5SummarizedExperiment(dir = here("processed-data",
                                                      "12_snRNA",
                                                      "sce_neuron_DimRed"))

print("Object loaded")
sce_neuron

#Additional sanity chekc
identical(rownames(colData(sce_neuron)),colnames(sce_neuron))
stopifnot(identical(rownames(colData(sce_neuron)),colnames(sce_neuron)))

sce_neuron

#######################################
#Begin clustering workflow.
print("Clustering")
Sys.time()
#build graph with k values of 10 and 50
print("Building shared nearest neighbor graph")
print("k=5")
snn_k_5  <- buildSNNGraph(sce_neuron,k = 5,use.dimred = "HARMONY",type  = "jaccard")
print("k=10")
snn_k_10 <- buildSNNGraph(sce_neuron, k = 10, use.dimred = "HARMONY",type="jaccard")
print("k=20")
snn_k_20 <- buildSNNGraph(sce_neuron, k = 20, use.dimred = "HARMONY",type="jaccard")
print("k=50")
snn_k_50 <- buildSNNGraph(sce_neuron, k = 50, use.dimred = "HARMONY",type="jaccard")  

#Louvain clustering
#k=5
set.seed(1234)
clust_5_1 <- igraph::cluster_louvain(snn_k_5,resolution=1)$membership
print("k=5, louvain res = 1")
table(clust_5_1)
sce_neuron$k_5_louvain_1 <- factor(clust_5_1)
k_5_1_tSNE <- plotReducedDim(object = sce_neuron,
                             dimred = "tSNE_HARMONY",
                             colour_by = "k_5_louvain_1",
                             text_by = "k_5_louvain_1") +
  ggtitle("k=5 louvain clustering\nresolution=1") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_5_1_tSNE, filename = here("plots","12_snRNA","Dim_Red","neuron_subset_k_5_louvain_res1_tSNE.png"))

#k=10
set.seed(1234)
clust_10_1 <- igraph::cluster_louvain(snn_k_10,resolution=1)$membership
print("k=10, louvain res = 1")
table(clust_10_1)
sce_neuron$k_10_louvain_1 <- factor(clust_10_1)
k_10_1_tSNE <- plotReducedDim(object = sce_neuron,
                              dimred = "tSNE_HARMONY",
                              colour_by = "k_10_louvain_1",
                              text_by = "k_10_louvain_1") +
  ggtitle("k=10 louvain clustering\nresolution=1") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_10_1_tSNE, filename = here("plots","12_snRNA","Dim_Red","neuron_subset_k_10_louvain_res1_tSNE.png"))

#Now k=20
set.seed(1234)
clust_20_1 <- igraph::cluster_louvain(snn_k_20,resolution=1)$membership
print("k=20, louvain res = 1")
table(clust_20_1)
sce_neuron$k_20_louvain_1 <- factor(clust_20_1)
k_20_1_tSNE <- plotReducedDim(object = sce_neuron,
                              dimred = "tSNE_HARMONY",
                              colour_by = "k_20_louvain_1",
                              text_by = "k_20_louvain_1") +
  ggtitle("k=20 louvain clustering\nresolution=1") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_20_1_tSNE, filename = here("plots","12_snRNA","Dim_Red","neuron_subset_k_20_louvain_res1_tSNE.png"))

#Now k=50
set.seed(1234)
clust_50_1 <- igraph::cluster_louvain(snn_k_50,resolution=1)$membership
print("k=50, louvain res = 1")
table(clust_50_1)
sce_neuron$k_50_louvain_1 <- factor(clust_50_1)
k_50_1_tSNE <- plotReducedDim(object = sce_neuron,
                              dimred = "tSNE_HARMONY",
                              colour_by = "k_50_louvain_1",
                              text_by = "k_50_louvain_1") +
  ggtitle("k=50 louvain clustering\nresolution=1") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_50_1_tSNE, filename = here("plots","12_snRNA","Dim_Red","neuron_subset_k_50_louvain_res1_tSNE.png"))

#Plot tSNE by sample
for(i in unique(sce_neuron$Brain_ID)){
  print(i)
  tSNE_BrainID <- plotReducedDim(object = sce_neuron[,sce_neuron$Brain_ID == i],
                                 dimred = "tSNE_HARMONY",
                                 colour_by = "Brain_ID") 
  ggsave(plot = tSNE_BrainID, filename = here("plots","12_snRNA",
                                              "Dim_Red","Neuron_Subset_BrainID",
                                               paste0(i,".png")))
}

print("Clustering complete")
#######################################
#Plot expression of known marker genes to be able to identify appropriate clusters
genes <- c("DRD1","DRD2","DRD3", #DA receptors
	   "RXFP1","CPNE4","CHST9","OPRM1","GRM8",#D1 islands
           "EBF1","HTR4", #D1 subsets from rat
	   "KIT","SST","CHAT","VIP", #Interneurons
           "RARB","CRYM", #striatal markers
	   "SLC17A7","ELAVL2","GAD1","GAD2","GLP1R") #misc

for(i in genes){
  print(i)
  x <- plotReducedDim(object = sce_neuron,
                      dimred = "tSNE_HARMONY",
                      colour_by = i,
                      swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","red")) +
  ggtitle(i) +
  theme(plot.title = element_text(hjust = 0.5))
  ggsave(plot = x,
         filename = here("plots","12_snRNA","Expression","Neuron_Subset",
                         paste0(i,".png")))
}


#Save the object
print("Saving sce object")
saveHDF5SummarizedExperiment(sce_neuron,
                             here("processed-data","12_snRNA",
                                  "sce_neuron_clustered"),
                             replace = TRUE)
print("Object saved")


print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
