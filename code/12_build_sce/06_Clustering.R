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

#Comoute log counts to plot expression
sce <- batchelor::multiBatchNorm(sce, batch = sce$Sample)

sce
# class: SingleCellExperiment 
# dim: 36601 118358 
# metadata(1): Samples
# assays(3): counts binomial_deviance_residuals logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(118358): 10_AAACCCAAGAAGTCAT-1 10_AAACCCACACGAGGAT-1 ...
# 9_TTTGGTTGTCAGCGTC-1 9_TTTGTTGAGATTAGTG-1
# colData names(22): Sample Barcode ... doubletScore sizeFactor
# reducedDimNames(4): GLMPCA_approx tSNE mnn tSNE_mnn
# mainExpName: NULL
# altExpNames(0):

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


#Lower K= fewer, larger clusters
#Smaller= more, small clusters.
#jaccard + louvain is similar to seurat workflow. 
#Resolution=1 is the default
snn_k_10 <- buildSNNGraph(sce, k = 10, use.dimred = "mnn",type="jaccard")

#Louvain clustering
set.seed(1234)
clust_10 <- igraph::cluster_louvain(snn_k_10,resolution=1)$membership
table(clust_10)
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 8412  5785  7813  7963 16092  6043  3871 10545  6057   989  8314  1289  2111 
# 14    15    16    17    18    19    20    21    22    23    24    25    26 
# 4515  3096  3740  1840  3492   301   488  1604   179  2632  1547   145   742 
# 27    28    29    30    31    32    33    34    35    36 
# 224   544    85   679   194  2030   989   472   993  2543 

#Add cluster information to object
sce$k_10_louvain <- factor(clust_10)

k_10_tSNE <- plotReducedDim(object = sce,
                            dimred = "tSNE_mnn",
                            colour_by = "k_10_louvain",
                            text_by = "k_10_louvain") +
  ggtitle("k=10 louvain clustering") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = k_10_tSNE, filename = here("plots","12_snRNA","Dim_Red","k_10_louvain_tSNE.png"))

#Goal right now is to identify if any non-neuronal cells are present within NeuN sorted samples. 
#We expect that some non-neuronal cells may get sorted but if a significant proportion of a NeuN 
#sorted sample is glia, then there could be an issue with staining. Either way, it will inform QC
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

#Classify each population as neuron or non-neuron. 
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

# 1 - Non
# 2 - Non
# 3 - Non
# 4 - Non
# 5 - Neuron
# 6 - Non
# 7 - Non
# 8 - Neuron
# 9 - Neuron
# 10 - Non
# 11 - Neuron
# 12 - Non
# 13 - Neuron
# 14 - NEuron
# 15 - Neuron
# 16 - Neuron
# 17 - NEuron
# 18 - Non
# 19 - Non
# 20 - Non
# 21 - Neuron
# 22 - Non
# 23 - Neuron
# 24 - Neuron
# 25 - Non
# 26 - Neuron
# 27 - Neuron
# 28 - Neuron
# 29 - Non
# 30 - Neuron
# 31 - Non 
# 32 - Neuron
# 33 - Neuron
# 34 - Neuron
# 35 - Neuron
# 36 - Neuron
#Non - 1,2,3,4,6,7,10,12,18,19,20,22,25,29,31
#Neuron - 5,8,9,11,13,14,15,16,17,21,23,24,26,27,28,30,32,33,34,35,36
sce$Binary_Classification <- ifelse(sce$k_10_louvain %in% c(5,8,9,11,13,14,15,16,17,21,23,
                                                            24,26,27,28,30,32,33,34,35,36),
                                    "Neuron",
                                    "Non-Neuron")

table(sce$Binary_Classification)
# Neuron Non-Neuron 
# 71309      47049 

#Color tSNE by Binary Classification
tSNE_bin <- plotReducedDim(object = sce,
                            dimred = "tSNE_mnn",
                            colour_by = "Binary_Classification") 
ggsave(plot = tSNE_bin, filename = here("plots","12_snRNA","Dim_Red","Binary_Classification_tSNE.png"))

#Plot # of genes per cluster and binary classification
genes_per_cluster   <- plotColData(object = sce,
                                   x = "k_10_louvain",
                                   y = "detected",
                                   colour_by = "k_10_louvain") +
    theme(legend.position = "none")
ggsave(filename = here("plots","12_snRNA","genes_per_cluster_violin.png"))
genes_per_bin_class <- plotColData(object = sce,
                                   x = "Binary_Classification",
                                   y = "detected",
                                   colour_by = "Binary_Classification") +
    theme(legend.position = "none")
ggsave(filename = here("plots","12_snRNA","genes_per_binary_classification_violin.png"))

#MAke a dataframe consisting of sample and binary classification information.
#as.data.frame.matrix maintains table structure 
x <- as.data.frame.matrix(table(sce$Sample,sce$Binary_Classification))
x$Sample <- row.names(x)

#Calculate the percentage of each sample that is 
x$Pct_Neuron <- (x$Neuron/(x$Neuron+x$`Non-Neuron`))*100
x$Pct_Non_Neuron <- (x$`Non-Neuron`/(x$Neuron+x$`Non-Neuron`))*100

#Factorize sample so that it is neuron first and non-neuron second
x$Sample <- factor(x = x$Sample,
                   levels = c("1c_NAc_SVB","3c_NAc_SVB","5c_NAc_SVB","7c_NAc_SVB",
                              "9c_NAc_SVB","11c_NAc_SVB","13c_NAc_SVB","15c_Nac_SVB",
                              "17c_Nac_SVB","19c_Nac_SVB","2c_NAc_SVB","4c_NAc_SVB",
                              "6c_NAc_SVB","8c_NAc_SVB","10c_NAc_SVB","12c_NAc_SVB",
                              "14c_NAc_SVB","16c_Nac_SVB","18c_Nac_SVB","20c_Nac_SVB"))

x_melt <- reshape2::melt(x[,c("Sample","Pct_Neuron","Pct_Non_Neuron")])

Pct_cell_bar <- ggplot(data = x_melt,aes(x = Sample,y = value, fill = variable)) +
    geom_bar(position="stack", stat="identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_segment(aes(x = 1, y = 102.5, xend = 10, yend = 102.5)) +
    geom_segment(aes(x = 11, y = 102.5, xend = 20, yend = 102.5)) +
    annotate(geom = "text",x = 5,y = 104,label = "PI+ Only") +
    annotate(geom = "text",x = 15,y = 104, label = "PI+NeuN+ Only")

ggsave(plot = Pct_cell_bar,filename = here("plots","12_snRNA","Pct_CellType_Barplot.png"))

