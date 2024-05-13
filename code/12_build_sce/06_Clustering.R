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

#Additional sanity chekc
identical(rownames(colData(sce)),colnames(sce))
stopifnot(identical(rownames(colData(sce)),colnames(sce)))

sce

#######################################
#Begin clustering workflow.

#build graph with k values of 10 and 50
print("Graph with k=20")
Sys.time()
snn_k_20 <- buildSNNGraph(sce, k = 20, use.dimred = "HARMONY")

####Also run walktrap clustering 
#k=20
print("Running walktrap clustering with k = 20")
Sys.time()

wt_clusters_k_20 <- igraph::cluster_walktrap(snn_k_20)$membership

print("k=20 done")
Sys.time()

#How many nuclei in each cluster?
table(wt_clusters_k_20)

#Add cluster information to the object
sce$k_20_walktrap <- factor(wt_clusters_k_20)

#tSNE with cluster information
k_20_wt_tSNE <- plotReducedDim(object = sce,
                              dimred = "tSNE_HARMONY",
                              colour_by = "k_20_walktrap",
                              text_by = "k_20_walktrap") +
  ggtitle("k=20 walktrap clustering") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_20_wt_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_20_walktrap_tSNE.png"))


#UMAP with cluster information
k_20_wt_umap <- plotReducedDim(object = sce,
                              dimred = "umap_HARMONY",
                              colour_by = "k_20_walktrap",
                              text_by = "k_20_walktrap") +
  ggtitle("k=20 walktrap clustering") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_20_wt_umap, filename = here("plots","12_snRNA","Dim_Red","k_20_walktrap_umap.png"))

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

print("Calculating log-normalized counts")
Sys.time()
sce <- computeSumFactors(sce,cluster = sce$k_20_walktrap,min.mean = 0.1)
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
           "MOBP","MBP","OPALIN", #Oligodendrocytes
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

#umap
for(i in genes){
  print(i)
  x <- plotReducedDim(object = sce,
                      dimred = "umap_HARMONY",
                      colour_by = i,
                      swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","red")) +
  ggtitle(i) +
  theme(plot.title = element_text(hjust = 0.5))
  ggsave(plot = x,
         filename = here("plots","12_snRNA","Expression","Known_Marker_Genes","UMAP",
                         paste0(i,"_UMAP.png")))
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
options(width = 120)
sessioninfo::session_info()
