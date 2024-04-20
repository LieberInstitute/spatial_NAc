#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

library(SingleCellExperiment)
library(sessioninfo)
library(HDF5Array)
library(ggplot2)
library(scater)
library(dplyr)
library(scran)
library(here)


##############
#Load sce
sce <- loadHDF5SummarizedExperiment(dir = here("processed-data","12_snRNA",
                                               "sce_clustered"))

dim(sce)
#[1]  36601 103785

identical(rownames(colData(sce)),colnames(sce))
#[1] TRUE

sce
# class: SingleCellExperiment 
# dim: 36601 103785 
# metadata(1): Samples
# assays(3): counts binomial_pearson_residuals logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(103785): 1_AAACCCAAGACCAACG-1 1_AAACCCACAGTCAGCC-1 ...
# 20_TTTGTTGCAAGATGTA-1 20_TTTGTTGGTACGAAAT-1
# colData names(39): Sample Barcode ... k_50_louvain_pt25 sizeFactor
# reducedDimNames(4): GLMPCA_approx tSNE HARMONY tSNE_HARMONY
# mainExpName: NULL
# altExpNames(0):

#load get mean ratio output to tryt and figure out what cluster 17 is
load(here("processed-data","12_snRNA","broad_clustering_meanratio_preannotation.rds"),
     verbose = TRUE)

seventeen <- subset(gmr_pt25,subset=(cellType.target == 17))

#Broad Cluster Annotations
# 1 - Oligo
# 2 - DRD1_MSN_A
# 3 - Microglia
# 4 - Polydendrocyte
# 5 - DRD2_MSN_A
# 6 - Ependymal
# 7 - Astrocyte
# 8 - DRD1_MSN_B
# 9 - Endothelial
# 10 - PVALB_Inh
# 11 - DRD2_MSN_B
# 12 - DRD1_MSN_C
# 13 - VIP_CCK_Inh
# 14 - GLP1R_Inh
# 15 - GABA_Glut
# 16 - SST_Inh
# 17 - T-Cells4

plotReducedDim(object = sce,
               dimred = "tSNE_HARMONY",
               colour_by = "SLC17A7",
               swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","red")) +
    ggtitle("SLC17A7") +
    theme(plot.title = element_text(hjust = 0.5))


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

#Annotate the tSNE
tsne_broad <- plotReducedDim(object = sce,
                             dimred = "tSNE_HARMONY",
                             colour_by = "CellType.Broad",
                             text_by = "CellType.Broad") +
    scale_color_manual(values = cluster_cols)

ggsave(filename = here("plots","12_snRNA","Dim_Red","tSNE_CellType_Broad.png"),plot = tsne_broad)

#########
#Calculate cluster percentage by BrainID
brain_by_cluster <- as.data.frame.matrix(table(sce$Brain_ID,sce$CellType.Broad))

#Calculate percentages 
brain_by_cluster_pct <- sweep(brain_by_cluster,MARGIN = 2,colSums(brain_by_cluster),"/") * 100
brain_by_cluster_pct$brain <- rownames(brain_by_cluster_pct)

#Melt dataframe and plot
brain_by_cluster_pct_melt <- reshape2::melt(brain_by_cluster_pct)

#make some new brain colors
brain_cols <- Polychrome::createPalette(length(unique(sce$Brain_ID)),
                                        c("#D81B60", "#1E88E5","#004D40"))
names(brain_cols) <- unique(sce$Brain_ID)

brain_cluster_bar <- ggplot(data = brain_by_cluster_pct_melt,aes(x = variable,y = value,fill = brain)) +
    geom_bar(position = "stack",stat = "identity") +
    scale_fill_manual(values = brain_cols) +
    labs(x = "Cluster",
         y = "Percent") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = brain_cluster_bar,
       filename = here("plots","12_snRNA","Dim_Red","BrainID_representation_by_cluster.png"))

###########
#Some QC measures. 
#Number of detected genes. 
detected_vln <- plotColData(object = sce,
                            y = "detected",
                            x = "CellType.Broad",
                            colour_by = "CellType.Broad") +
    scale_color_manual(values = cluster_cols) +
    ggtitle("Number of Detected Genes") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) 
ggsave(plot = detected_vln,
       filename = here("plots","12_snRNA","Broad_CellType_detected.png"))

#Library size 
library_size_vln <- plotColData(object = sce,
                                y = "sum",
                                x = "CellType.Broad",
                                colour_by = "CellType.Broad") +
    scale_y_log10() +
    scale_color_manual(values = cluster_cols) +
    ggtitle("Library Size") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) 
ggsave(plot = library_size_vln,
       filename = here("plots","12_snRNA","Broad_CellType_library_size.png"))

#percent mito
percent_mito_vln <- plotColData(object = sce,
                                y = "subsets_Mito_percent",
                                x = "CellType.Broad",
                                colour_by = "CellType.Broad") +
    scale_color_manual(values = cluster_cols) +
    ggtitle("Percent mito") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) 
ggsave(plot = percent_mito_vln,
       filename = here("plots","12_snRNA","Broad_CellType_PercentMito.png"))

#doublet score
doublet_boxplot <- colData(sce) %>% 
    as.data.frame() %>%
    select(CellType.Broad,doubletScore) %>%
    ggplot(aes(x = CellType.Broad,y = doubletScore,fill = CellType.Broad)) +
    geom_boxplot() +
    scale_fill_manual(values = cluster_cols) +
    geom_hline(yintercept = 5,lty = 2,color = "red") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = doublet_boxplot,
       filename = here("plots","12_snRNA","Broad_CellType_doubletScore.png"))

#Number of cells barplot
Num_Nucs <- table(sce$CellType.Broad) %>% as.data.frame() %>% 
    ggplot(aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cluster_cols) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Broad Cell Type",
         y = "Number of Nuclei") 

ggsave(plot = Num_Nucs,
       filename = here("plots","12_snRNA","Broad_CellType_NumberofNuclei_Bargraph.png"))

# #Get the cell IDs of the IEG_Neurons to remove them. 
# IEG_cells <- rownames(colData(sce)[which(sce$CellType.Broad == "IEG_neurons"),])
# save(IEG_cells,file = here("processed-data","12_snRNA","IEG_neurons_to_remove.rds"))

#Make a colData column for neurons or glia
sce$Neuron_NonNeuron <- ifelse(sce$CellType.Broad %in% c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C",
                                                         "DRD2_MSN_A","DRD2_MSN_B","GLP1R_Inh",
                                                         "PVALB_Inh","VIP_CCK_Inh","SST_Inh",
                                                         "GABA_Glut"),
                               "Neuron",
                               "NonNeuron")

table(sce$Neuron_NonNeuron)
# Neuron NonNeuron 
# 66042     37743


Neuron_NonNeuron_bar <- table(sce$Neuron_NonNeuron) %>% as.data.frame() %>% 
    ggplot(aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat = "identity") +
    labs(x = "Cell Class",
         y = "Number of Nuclei") +
    theme_bw() +
    theme(legend.position = "none") 
ggsave(plot = Neuron_NonNeuron_bar,
       filename = here("plots","12_snRNA","NeuronNonNeuron_Bargraph.png"))


#Load in the object with QC metrics to make some additional plots showing which cells are removed. 
load(here("processed-data","12_snRNA","sce_emptyDrops_removed_withQC.rds"),verbose = TRUE)

sce
# class: SingleCellExperiment 
# dim: 36601 120449 
# metadata(1): Samples
# assays(1): counts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(6): source type ... gene_name gene_type
# colnames(120449): 1_AAACCCAAGACCAACG-1 1_AAACCCACAGTCAGCC-1 ...
# 20_TTTGTTGCAAGATGTA-1 20_TTTGTTGGTACGAAAT-1
# colData names(30): Sample Barcode ... IEG_cells discard
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

#library size
lib_size_qc <- plotColData(object = sce,y = "sum",x = "Sample",colour_by = "discard_sum_detected") +
    scale_y_continuous(trans = "log10") +
    ggtitle("Library Size") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) 
ggsave(plot = lib_size_qc,filename = here("plots","12_snRNA","nuclei_QC","lib_size_QC.png"))

#detected features
genes_qc <- plotColData(object = sce,y = "detected",x = "Sample",colour_by = "discard_sum_detected") +
    scale_y_continuous(trans = "log10") +
    ggtitle("Detected Features") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) 
ggsave(plot = genes_qc,filename = here("plots","12_snRNA","nuclei_QC","detected_features_QC.png"))
