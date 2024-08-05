# cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
# module load r_nac

#Load libraries
library(SingleCellExperiment)
library(sessioninfo)
library(batchelor)
library(rafalib)
library(ggplot2)
library(scater)
library(dplyr)
library(scran)
library(here)

#load sce

print("Loading Object")
Sys.time()
sce <- readRDS(here("processed-data","12_snRNA","sce_DimRed.Rds"))
print("Object loaded")
Sys.time()

sce

#Additional sanity check
identical(rownames(colData(sce)),colnames(sce))
stopifnot(identical(rownames(colData(sce)),colnames(sce)))

sce

#######################################
#Begin clustering workflow.
##############K=5####################
#build graph with k value of 5
print("Graph with k=5")
Sys.time()
set.seed(5)
snn_k_5 <- buildSNNGraph(sce, k = 5, use.dimred = "HARMONY")

#louvain clustering with resolution of 0.5
print("Running louvain k=5, res 0.5")
Sys.time()
set.seed(5)
louvain_clusters_k_5_pt5 <- igraph::cluster_louvain(snn_k_5,resolution = 0.5)$membership
table(louvain_clusters_k_5_pt5)
Sys.time()

#louvain clustering with resolution of 1
print("Running louvain k=5, res 1")
Sys.time()
set.seed(5)
louvain_clusters_k_5_1 <- igraph::cluster_louvain(snn_k_5,resolution = 1)$membership
table(louvain_clusters_k_5_1)
Sys.time()

#Add cluster information to the object
sce$k_5_louvain_pt5 <- factor(louvain_clusters_k_5_pt5)
sce$k_5_louvain_1 <- factor(louvain_clusters_k_5_1)

#tSNE with cluster information
k_5_pt5_tSNE <- plotReducedDim(object = sce,
                                dimred = "tSNE_HARMONY",
                                colour_by = "k_5_louvain_pt5",
                                text_by = "k_5_louvain_pt5") +
  ggtitle("k=5 louvain clustering, resolution = 0.5") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_5_pt5_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_5_louvain_pt5_tSNE.png"))

#tSNE with cluster information
k_5_1_tSNE <- plotReducedDim(object = sce,
                              dimred = "tSNE_HARMONY",
                              colour_by = "k_5_louvain_1",
                              text_by = "k_5_louvain_1") +
  ggtitle("k=5 louvain clustering, resolution = 1") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_5_1_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_5_louvain_1_tSNE.png"))



##############K=10####################
#build graph with k value of 10
print("Graph with k=10")
Sys.time()
set.seed(10)
snn_k_10 <- buildSNNGraph(sce, k = 10, use.dimred = "HARMONY")

#louvain clustering with resolution of 0.5
print("Running louvain k=10, res 0.5")
Sys.time()
set.seed(10)
louvain_clusters_k_10_pt5 <- igraph::cluster_louvain(snn_k_10,resolution = 0.5)$membership
table(louvain_clusters_k_10_pt5)
Sys.time()

#louvain clustering with resolution of 1
print("Running louvain k=10, res 1")
Sys.time()
set.seed(10)
louvain_clusters_k_10_1 <- igraph::cluster_louvain(snn_k_10,resolution = 1)$membership
table(louvain_clusters_k_10_1)
Sys.time()

#Add cluster information to the object
sce$k_10_louvain_pt5 <- factor(louvain_clusters_k_10_pt5)
sce$k_10_louvain_1 <- factor(louvain_clusters_k_10_1)

#tSNE with cluster information
k_10_pt5_tSNE <- plotReducedDim(object = sce,
                                dimred = "tSNE_HARMONY",
                                colour_by = "k_10_louvain_pt5",
                                text_by = "k_10_louvain_pt5") +
  ggtitle("k=10 louvain clustering, resolution = 0.5") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_10_pt5_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_10_louvain_pt5_tSNE.png"))

#tSNE with cluster information
k_10_1_tSNE <- plotReducedDim(object = sce,
                              dimred = "tSNE_HARMONY",
                              colour_by = "k_10_louvain_1",
                              text_by = "k_10_louvain_1") +
  ggtitle("k=10 louvain clustering, resolution = 1") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_10_1_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_10_louvain_1_tSNE.png"))


##############K=20####################
#build graph with k value of 20
print("Graph with k=20")
Sys.time()
set.seed(20)
snn_k_20 <- buildSNNGraph(sce, k = 20, use.dimred = "HARMONY")

#louvain clustering with resolution of 0.5
print("Running louvain k=20, res 0.5")
Sys.time()
set.seed(20)
louvain_clusters_k_20_pt5 <- igraph::cluster_louvain(snn_k_20,resolution = 0.5)$membership
table(louvain_clusters_k_20_pt5)
Sys.time()

#louvain clustering with resolution of 1
print("Running louvain k=20, res 1")
Sys.time()
set.seed(20)
louvain_clusters_k_20_1 <- igraph::cluster_louvain(snn_k_20,resolution = 1)$membership
table(louvain_clusters_k_20_1)
Sys.time()

#Add cluster information to the object
sce$k_20_louvain_pt5 <- factor(louvain_clusters_k_20_pt5)
sce$k_20_louvain_1 <- factor(louvain_clusters_k_20_1)

#tSNE with cluster information
k_20_pt5_tSNE <- plotReducedDim(object = sce,
                                dimred = "tSNE_HARMONY",
                                colour_by = "k_20_louvain_pt5",
                                text_by = "k_20_louvain_pt5") +
  ggtitle("k=20 louvain clustering, resolution = 0.5") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_20_pt5_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_20_louvain_pt5_tSNE.png"))

#tSNE with cluster information
k_20_1_tSNE <- plotReducedDim(object = sce,
                                dimred = "tSNE_HARMONY",
                                colour_by = "k_20_louvain_1",
                                text_by = "k_20_louvain_1") +
  ggtitle("k=20 louvain clustering, resolution = 1") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_20_1_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_20_louvain_1_tSNE.png"))


##############K=50####################
#build graph with k value of 50
print("Graph with k=50")
Sys.time()
set.seed(50)
snn_k_50 <- buildSNNGraph(sce, k = 50, use.dimred = "HARMONY")

#louvain clustering with resolution of 0.5
print("Running louvain k=50, res 0.5")
Sys.time()
set.seed(50)
louvain_clusters_k_50_pt5 <- igraph::cluster_louvain(snn_k_50,resolution = 0.5)$membership
table(louvain_clusters_k_50_pt5)
Sys.time()

#louvain clustering with resolution of 1
print("Running louvain k=50, res 1")
Sys.time()
set.seed(50)
louvain_clusters_k_50_1 <- igraph::cluster_louvain(snn_k_50,resolution = 1)$membership
table(louvain_clusters_k_50_1)
Sys.time()

#Add cluster information to the object
sce$k_50_louvain_pt5 <- factor(louvain_clusters_k_50_pt5)
sce$k_50_louvain_1 <- factor(louvain_clusters_k_50_1)

#tSNE with cluster information
k_50_pt5_tSNE <- plotReducedDim(object = sce,
                                dimred = "tSNE_HARMONY",
                                colour_by = "k_50_louvain_pt5",
                                text_by = "k_50_louvain_pt5") +
  ggtitle("k=50 louvain clustering, resolution = 0.5") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_50_pt5_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_50_louvain_pt5_tSNE.png"))

#tSNE with cluster information
k_50_1_tSNE <- plotReducedDim(object = sce,
                              dimred = "tSNE_HARMONY",
                              colour_by = "k_50_louvain_1",
                              text_by = "k_50_louvain_1") +
  ggtitle("k=50 louvain clustering, resolution = 1") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_50_1_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_50_louvain_1_tSNE.png"))

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
Sys.time()
#######################################

print("Calculating log-normalized counts") #Going to jsut use the k=10, louvain 1 clustering --> Will use this clustering to generate final cell types
Sys.time()
sce <- computeSumFactors(sce,cluster = sce$k_10_louvain_1,min.mean = 0.1)
sce <- logNormCounts(sce)
Sys.time()

#Plot expression of known marker genes to be able to identify appropriate clusters
genes <- c("SNAP25","SYT1","RBFOX3", #PAN-NEURON markers
           "GAD1","GAD2","SLC32A1","ELAVL2", #GABA markers
           "GPR88","NEXN", #BROAD STRIATAL MARKER FROM DROPVIZ
           "OTOF","CASZ1","NRIP3","PCDH8", #BNST markers from dropviz
           "PPP1R1B","BCL11B","ISL1","FOXP2","RELN","RARB",#MSN markers
           "DRD1","PDYN","EBF1","RXFP1","TAC1","CRHR2","RXFP1","CPNE4",#D1 markers
           "DRD2","PENK","ADORA2A",#D2 markers,
           "GRM8","CHST9","OPRM1", #D1-ISLANDS
           "DRD3","KCNT2","NTN1","TRPM3","PLD5", #ICJ
           "KIT","CHAT","SST", #INTERNEURON MARKERS
           "SLC18A3","SLC5A7", #Additional CHAT markers 
           "MOBP","MBP","OPALIN", #Oligodendrocytes
           "SLC17A6","SLC17A7", #excitatory
           "PDGFRA", #POLYDENDROCYTES
           "AQP4","GFAP", #Astrocytes
           "CD74", #Microglia
           "TAC3","GLP1R") #Misc. markers

#tSNE
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
         filename = here("plots","12_snRNA","Expression","Known_Marker_Genes","tSNE",
                         paste0(i,"_tSNE.png")))
}


#Violin plots by each clustering 
for(i in genes){
  for(l in c("k_5_louvain_pt5","k_5_louvain_1",
             "k_10_louvain_pt5","k_10_louvain_1",
             "k_20_louvain_pt5","k_20_louvain_1",
             "k_50_louvain_pt5","k_50_louvain_1")){
    y <- plotExpression(sce,
                        x = l, #x-axis by clustering 
                        features = i, #feature is one of genes listed above
                        swap_rownames = "gene_name") +
      theme(axis.text.x = element_text(angle=90,hjust = 1),
            legend.position = "none") +
      stat_summary(fun = median, 
                   fun.min = median, 
                   fun.max = median,
                   geom = "crossbar", 
                   width = 0.3) 
    ggsave(plot = y,
           filename = here("plots","12_snRNA","Expression",
                           "Known_Marker_Genes","Violin_Plots",l,
                           paste0(i,l,"_Violin.png")))
  }
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
#######################################

#Save the object
print("Saving sce object")
Sys.time()
saveRDS(sce,file = here("processed-data","12_snRNA","sce_clustered_log.Rds"))

print("Object saved")

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 150)
sessioninfo::session_info()
