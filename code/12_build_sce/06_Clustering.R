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

sce

#######################################
#Begin clustering workflow.
print("Clustering")
Sys.time()
#build graph with k values of 10 and 50
snn_k_10 <- buildSNNGraph(sce, k = 10, use.dimred = "HARMONY",type="jaccard")
snn_k_50 <- buildSNNGraph(sce, k = 50, use.dimred = "HARMONY",type="jaccard")  

#Louvain clustering
#Try different resolution values as well. Lower resolution values result in larger, fewer clusters
#First K=10
set.seed(1234)
#1
clust_10_1 <- igraph::cluster_louvain(snn_k_10,resolution=1)$membership
print("k=10, louvain res = 1")
table(clust_10_1)
sce$k_10_louvain_1 <- factor(clust_10_1)
k_10_1_tSNE <- plotReducedDim(object = sce,
                              dimred = "tSNE_HARMONY",
                              colour_by = "k_10_louvain_1",
                              text_by = "k_10_louvain_1") +
  ggtitle("k=10 louvain clustering\nresolution=1") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_10_1_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_10_louvain_res1_tSNE.png"))

#0.75
clust_10.75 <- igraph::cluster_louvain(snn_k_10,resolution=0.75)$membership
print("k=10, louvain res = .75")
table(clust_10.75)
sce$k_10_louvain_pt75 <- factor(clust_10.75)
k_10_pt75_tSNE <- plotReducedDim(object = sce,
                            dimred = "tSNE_HARMONY",
                            colour_by = "k_10_louvain_pt75",
                            text_by = "k_10_louvain_pt75") +
  ggtitle("k=10 louvain clustering\nresolution=0.75") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_10_pt75_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_10_louvain_respt75_tSNE.png"))

#0.5
clust_10.5 <- igraph::cluster_louvain(snn_k_10,resolution=0.5)$membership
print("k=10, louvain res = .5")
table(clust_10.5)
sce$k_10_louvain_pt5 <- factor(clust_10.5)
k_10_pt5_tSNE <- plotReducedDim(object = sce,
                                dimred = "tSNE_HARMONY",
                                colour_by = "k_10_louvain_pt5",
                                text_by = "k_10_louvain_pt5") +
  ggtitle("k=10 louvain clustering\nresolution=0.5") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_10_pt5_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_10_louvain_respt5_tSNE.png"))

#0.25
clust_10.25 <- igraph::cluster_louvain(snn_k_10,resolution=0.25)$membership
sce$k_10_louvain_pt25 <- factor(clust_10.25)
table(clust_10.25)
k_10_pt25_tSNE <- plotReducedDim(object = sce,
                                 dimred = "tSNE_HARMONY",
                                 colour_by = "k_10_louvain_pt25",
                                 text_by = "k_10_louvain_pt25") +
  ggtitle("k=10 louvain clustering\nresolution=0.25") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_10_pt25_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_10_louvain_respt25_tSNE.png"))

#Now k=50
set.seed(1234)
#1
clust_50_1 <- igraph::cluster_louvain(snn_k_50,resolution=1)$membership
print("k=50, louvain res = 1")
table(clust_50_1)
sce$k_50_louvain_1 <- factor(clust_50_1)
k_50_1_tSNE <- plotReducedDim(object = sce,
                              dimred = "tSNE_HARMONY",
                              colour_by = "k_50_louvain_1",
                              text_by = "k_50_louvain_1") +
  ggtitle("k=50 louvain clustering\nresolution=1") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_50_1_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_50_louvain_res1_tSNE.png"))

#0.75
clust_50.75 <- igraph::cluster_louvain(snn_k_50,resolution=0.75)$membership
print("k=50, louvain res = .75")
table(clust_50.75)
sce$k_50_louvain_pt75 <- factor(clust_50.75)
k_50_pt75_tSNE <- plotReducedDim(object = sce,
                            dimred = "tSNE_HARMONY",
                            colour_by = "k_50_louvain_pt75",
                            text_by = "k_50_louvain_pt75") +
  ggtitle("k=50 louvain clustering\nresolution=0.75") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_50_pt75_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_50_louvain_respt75_tSNE.png"))

#0.5
clust_50.5 <- igraph::cluster_louvain(snn_k_50,resolution=0.5)$membership
print("k=50, louvain res = .5")
table(clust_50.5)
sce$k_50_louvain_pt5 <- factor(clust_50.5)
k_50_pt5_tSNE <- plotReducedDim(object = sce,
                                dimred = "tSNE_HARMONY",
                                colour_by = "k_50_louvain_pt5",
                                text_by = "k_50_louvain_pt5") +
  ggtitle("k=50 louvain clustering\nresolution=0.5") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_50_pt5_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_50_louvain_respt5_tSNE.png"))

#0.25
clust_50.25 <- igraph::cluster_louvain(snn_k_50,resolution=0.25)$membership
sce$k_50_louvain_pt25 <- factor(clust_50.25)
table(clust_50.25)
k_50_pt25_tSNE <- plotReducedDim(object = sce,
                                 dimred = "tSNE_HARMONY",
                                 colour_by = "k_50_louvain_pt25",
                                 text_by = "k_50_louvain_pt25") +
  ggtitle("k=50 louvain clustering\nresolution=0.25") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_50_pt25_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_50_louvain_respt25_tSNE.png"))


#Plot tSNE by sample
for(i in unique(sce$Sample)){
  print(i)
  tSNE_sample <- plotReducedDim(object = sce[,sce$Sample == i],
                                dimred = "tSNE_HARMONY",
                                colour_by = "Sample") +
    ggtitle(unique(colData(sce[,sce$Sample == i])[,"Sort"])) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(plot = tSNE_sample, filename = here("plots","12_snRNA",
                                             "Dim_Red","tSNE_Sample_Specific",
                                             paste0(i,".png")))
}

print("Clustering complete")
#######################################

#Calcualte log counts and check expression of known marker genes to aid in cluster annotation
#use k=10, 0.25 resolution for log-normalization. Will have to recalculate if choose different clustering. 
print("Calculating log-normalized counts")
Sys.time()
sce <- computeSumFactors(sce,cluster = sce$k_10_louvain_pt25,min.mean = 0.1)
sce <- logNormCounts(sce)

#Plot expression of known marker genes to be able to identify appropriate clusters
genes <- c("SNAP25","SYT1","RBFOX3", #PAN-NEURON markers
           "GAD1","GAD2","SLC32A1","ELAVL2", #GABA markers
           "GPR88","NEXN", #BROAD STRIATAL MARKER FROM DROPVIZ
           "OTOF","CASZ1","NRIP3","PCDH8", #BNST markers from dropviz
           "PPP1R1B","BCL11B","ISL1","FOXP2","RELN","RARB",#MSN markers
           "DRD1","PDYN","EBF1","RXFP1","TAC1","CRHR2","RXFP1","CPNE4",#D1 markers
           "DRD2","PENK","ADORA2A",#D2 markers,
           "GRM8","CHST9","OPRM1", #D1-ISLANDS
           "KIT","CHAT","SST", #INTERNEURON MARKERS
           "MOBP","MBP","OPALIN", #Oligodendrocytes
           "PDGFRA", #POLYDENDROCYTES
           "AQP4","GFAP", #Astrocytes
           "CD74", #Microglia
           "TAC3","GLP1R") #Misc. markers

for(i in genes){
  print(i)
  x <- plotReducedDim(object = sce,
                      dimred = "tSNE_HARMONY",
                      colour_by = i,
                      swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","red")) +
  ggtitle(i) +
  theme(plot.title = element_text(hjust = 0.5))
  ggsave(plot = x,
         filename = here("plots","12_snRNA","Expression","Known_Marker_Genes",
                         paste0(i,".png")))
}

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

#Save the object
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
