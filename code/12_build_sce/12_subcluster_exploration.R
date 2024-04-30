#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

#load libraries
library(SingleCellExperiment)
library(sessioninfo)
library(HDF5Array)
library(ggplot2)
library(scater)
library(dplyr)
library(scran)
library(here)

#load sce
sce_neuron <- loadHDF5SummarizedExperiment(dir = here("processed-data",
                                                      "12_snRNA",
                                                      "sce_neuron_clustered"))

sce_neuron
# class: SingleCellExperiment 
# dim: 36601 66042 
# metadata(1): Samples
# assays(3): counts binomial_pearson_residuals logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(66042): 1_AAACCCATCACCCTCA-1 1_AAACGAACATAGACTC-1 ...
# 20_TTTGTTGCAAGATGTA-1 20_TTTGTTGGTACGAAAT-1
# colData names(43): Sample Barcode ... k_5_louvain_1 k_20_louvain_1
# reducedDimNames(4): GLMPCA tSNE HARMONY tSNE_HARMONY
# mainExpNameNULL
# altExpNames(0):


#Calculate cluster percentage by BrainID
brain_by_cluster <- as.data.frame.matrix(table(sce_neuron$Brain_ID,sce_neuron$k_50_louvain_1))

#Calculate percentages 
brain_by_cluster_pct <- sweep(brain_by_cluster,MARGIN = 2,colSums(brain_by_cluster),"/") * 100
brain_by_cluster_pct$brain <- rownames(brain_by_cluster_pct)

#Melt dataframe and plot
brain_by_cluster_pct_melt <- reshape2::melt(brain_by_cluster_pct)

#make some new brain colors
brain_cols <- Polychrome::createPalette(length(unique(sce_neuron$Brain_ID)),
                                        c("#D81B60", "#1E88E5","#004D40"))
names(brain_cols) <- unique(sce_neuron$Brain_ID)

brain_cluster_bar <- ggplot(data = brain_by_cluster_pct_melt,aes(x = variable,y = value,fill = brain)) +
    geom_bar(position = "stack",stat = "identity") +
    scale_fill_manual(values = brain_cols) +
    labs(x = "Cluster",
         y = "Percent") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = brain_cluster_bar,filename = here("plots","12_snRNA",
                                                "Dim_Red","Neuron_Subset_BrainID",
                                                "k_50_louvain_by_brainID.png"))

#Evaluate expression of several genes within the subclusters. 
subset_marker_gene_vln <- plotExpression(object = sce_neuron,
                                         features = c("DRD1","DRD2","DRD3",
                                                      "RXFP1","CPNE4","OPRM1",
                                                      "PVALB","KIT","SST",
                                                      "CHAT","GAD1","SLC17A7"),
                                         colour_by = "k_50_louvain_1",
                                         x = "k_50_louvain_1",
                                         ncol = 3,
                                         swap_rownames = "gene_name") +
    stat_summary(fun = median, geom = "crossbar", width = 0.3)

ggsave(plot = subset_marker_gene_vln,filename = here("plots","12_snRNA",
                                                "Expression","Neuron_Subset",
                                              "subset_marker_gene_vln.png"),
       height = 10.5,width = 20)

celltypebroad_marker_gene_vln <- plotExpression(object = sce_neuron,
                                                features = c("DRD1","DRD2","DRD3",
                                                             "RXFP1","CPNE4","OPRM1",
                                                             "PVALB","KIT","SST",
                                                             "CHAT","GAD1","SLC17A7"),
                                                colour_by = "CellType.Broad",
                                                x = "CellType.Broad",
                                                ncol = 3,
                                                swap_rownames = "gene_name") +
    theme(axis.text.x = element_text(angle = 45,hjust =1)) +
    stat_summary(fun = median, geom = "crossbar", width = 0.3)

ggsave(plot = celltypebroad_marker_gene_vln,filename = here("plots","12_snRNA",
                                                            "Expression",
                                                            "celltype_broad_marker_gene_vln.png"),
       height = 10.5,width = 20)

#Remove levels of CellType.Broad that map to non-neuronal populations.
sce_neuron$CellType.Broad <- factor(x = sce_neuron$CellType.Broad,
                                    levels = levels(sce_neuron$CellType.Broad)[1:10])

#see the subclusters that make up the broad cell types
sub_by_cluster <- as.data.frame.matrix(table(sce_neuron$CellType.Broad,sce_neuron$k_50_louvain_1))

#Calculate cluster percentage by BrainID
sub_by_cluster_pct <- sweep(sub_by_cluster,MARGIN = 2,colSums(sub_by_cluster),"/") * 100
sub_by_cluster_pct$subcluster_k_50 <- rownames(sub_by_cluster_pct)

#Melt dataframe and plot
sub_by_cluster_pct_melt <- reshape2::melt(sub_by_cluster_pct)

#make some new brain colors
sub_cols <- Polychrome::createPalette(length(unique(sce_neuron$CellType.Broad)),
                                        c("#D81B60", "#1E88E5","#004D40"))
names(sub_cols) <- unique(sce_neuron$CellType.Broad)

sub_cluster_bar <- ggplot(data = sub_by_cluster_pct_melt,aes(x = variable,y = value,fill = suluster_k_50)) +
    geom_bar(position = "stack",stat = "identity") +
    scale_fill_manual(values = sub_cols) +
    labs(x = "Subcluster (k=50)",
         y = "Percent") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = sub_cluster_bar,filename = here("plots","12_snRNA",
                                                "Dim_Red",
                                                "BroadCellType_in_Subcluster.png"))
