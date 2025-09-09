#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

library(SingleCellExperiment)
library(ComplexHeatmap)
library(sessioninfo)
library(ggplot2)
library(scater)
library(dplyr)
library(scran)
library(here)

#load the sce object
sce <- readRDS(here("processed-data","12_snRNA","sce_CellType_noresiduals.Rds"))

sce
# class: SingleCellExperiment 
# dim: 36601 103785 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(103785): 1_AAACCCAAGACCAACG-1 1_AAACCCACAGTCAGCC-1 ...
# 20_TTTGTTGCAAGATGTA-1 20_TTTGTTGGTACGAAAT-1
# colData names(41): Sample Barcode ... sizeFactor CellType.Final
# reducedDimNames(4): GLMPCA_approx tSNE HARMONY tSNE_HARMONY
# mainExpName: NULL
# altExpNames(0):

#Load cluster colors
load(here("processed-data","12_snRNA","070924_21colors_celltypeFinal.rda"),verbose = TRUE)
# Loading objects:
#   cluster_cols

#Save the tSNE as a pdf 
tSNE_HARMONY <- plotReducedDim(object = sce,
                               dimred = "tSNE_HARMONY",
                               colour_by = "CellType.Final",
                               text_by = "CellType.Final") +
  scale_color_manual(values = cluster_cols) +
  theme(legend.position = "none")

ggsave(plot = tSNE_HARMONY,filename = here("plots","12_snRNA","tSNE_CellType_Final.pdf"),
       height = 10, width = 10)


#color the tSNE by sort
tSNE_Sort <- plotReducedDim(object = sce,
                               dimred = "tSNE_HARMONY",
                               colour_by = "Sort")

ggsave(plot = tSNE_Sort,filename = here("plots","12_snRNA","tSNE_ColoredbySort_Fig1.pdf"),
       height = 10, width = 10)


plotExpression(object = sce,features = c("KIT","CHAT","SST","VIP","CORT","PNOC"),x = "CellType.Final",
               colour_by = "CellType.Final",swap_rownames = "gene_name",ncol = 2) +
  scale_color_manual(values = cluster_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)

###Feature plots for SNAP25, GAD1, DRD1, and DRD2
#SYT1
syt1_featureplot <- plotReducedDim(sce,
                                     dimred = "tSNE_HARMONY",
                                     colour_by = "SYT1",
                                     swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","red"))

ggsave(plot = syt1_featureplot,
       filename = here("plots","12_snRNA","Expression","SYT1_Featureplot_Figure1.pdf"),
       height = 8,width = 8)

#SNAP25
snap25_featureplot <- plotReducedDim(sce,
                                   dimred = "tSNE_HARMONY",
                                   colour_by = "SNAP25",
                                   swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","red"))

ggsave(plot = snap25_featureplot,
       filename = here("plots","12_snRNA","Expression","SNAP25_Featureplot_Figure1.pdf"),
       height = 8,width = 8)


#RBFOX3
rbfox3_featureplot <- plotReducedDim(sce,
                                   dimred = "tSNE_HARMONY",
                                   colour_by = "RBFOX3",
                                   swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","red"))

ggsave(plot = rbfox3_featureplot,
       filename = here("plots","12_snRNA","Expression","RBFOX3_Featureplot_Figure1.pdf"),
       height = 8,width = 8)

#GAD1
gad1_featureplot <- plotReducedDim(sce,
                                     dimred = "tSNE_HARMONY",
                                     colour_by = "GAD1",
                                     swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","red"))

ggsave(plot = gad1_featureplot,
       filename = here("plots","12_snRNA","Expression","GAD1_Featureplot_Figure1.pdf"),
       height = 8,width = 8)

#DRD1
D1_featureplot <- plotReducedDim(sce,
                                 dimred = "tSNE_HARMONY",
                                 colour_by = "DRD1",
                                 swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","red"))

ggsave(plot = D1_featureplot,
       filename = here("plots","12_snRNA","Expression","DRD1_Featureplot_Figure1.pdf"),
       height = 8,width = 8)

#DRD2
D2_featureplot <- plotReducedDim(sce,
                                 dimred = "tSNE_HARMONY",
                                 colour_by = "DRD2",
                                 swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","red"))

ggsave(plot = D2_featureplot,
       filename = here("plots","12_snRNA","Expression","DRD2_Featureplot_Figure1.pdf"),
       height = 8,width = 8)

#Remove Neuron_Ambig from further analysis. 
sce <- sce[,sce$CellType.Final != "Neuron_Ambig"]

sce
# class: SingleCellExperiment 
# dim: 36601 103339 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(103339): 1_AAACCCAAGACCAACG-1 1_AAACCCACAGTCAGCC-1 ...
# 20_TTTGTTGCAAGATGTA-1 20_TTTGTTGGTACGAAAT-1
# colData names(41): Sample Barcode ... sizeFactor CellType.Final
# reducedDimNames(4): GLMPCA_approx tSNE HARMONY tSNE_HARMONY
# mainExpName: NULL
# altExpNames(0):

#Make order for CellType.Final 
sce$CellType.Final <- factor(x = sce$CellType.Final,
                             levels = c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                                        "DRD2_MSN_A","DRD2_MSN_B",
                                        "Inh_A","Inh_B","Inh_C","Inh_D","Inh_E","Inh_F",
                                        "Excitatory",
                                        "Astrocyte_A","Astrocyte_B","Ependymal",
                                        "Oligo","OPC",
                                        "Microglia",
                                        "Endothelial"))


####Remake BrainID by cluster with new order. 
#Calculate cluster percentage by BrainID
brain_by_cluster <- as.data.frame.matrix(table(sce$Brain_ID,sce$CellType.Final))

#Calculate percentages 
brain_by_cluster_pct <- sweep(brain_by_cluster,MARGIN = 2,STATS = colSums(brain_by_cluster),FUN = "/") * 100
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
  labs(x = "CellType",
       y = "Percent") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = brain_cluster_bar,
       filename = here("plots","12_snRNA","Dim_Red","BrainID_by_CellType_Final_Ordered.pdf"),
       height = 12,
       width = 12)

############
#Remake the above plot with donor on the x-axis. 
cluster_by_brain <- as.data.frame.matrix(table(sce$CellType.Final,sce$Brain_ID))

#Calculate percentages 
cluster_by_brain_pct <- sweep(cluster_by_brain,MARGIN = 2,STATS = colSums(cluster_by_brain),FUN = "/") * 100
cluster_by_brain_pct$CellType <- rownames(cluster_by_brain_pct)

#Melt dataframe and plot
cluster_by_brain_pct_melt <- reshape2::melt(cluster_by_brain_pct)
#Using CellType as id variables

cluster_cols
# Oligo  DRD1_MSN_A  DRD2_MSN_A         OPC   Microglia   Ependymal 
# "#4F4753"   "#ECA31C"   "#58B6ED"   "#0D9F72"   "#F2E642"   "#0077B9" 
# Astrocyte_A  DRD1_MSN_B Endothelial       Inh_A  DRD2_MSN_B Astrocyte_B 
# "#D95F00"   "#D079AA"   "#D00DFF"   "#35FB00"   "#F80091"   "#FF0016" 
# DRD1_MSN_C Neuro_Ambig  DRD1_MSN_D       Inh_B       Inh_C       Inh_D 
# "#2A4BF9"      "grey"   "#FB3DD9"   "#7A0096"   "#854222"   "#A7F281" 
# Inh_E  Excitatory       Inh_F 
# "#0DFBFA"   "#5C6300"     "black" 

#Remove Neuron_Ambig
cluster_cols <- cluster_cols[-14]

cluster_brain_bar <- ggplot(data = cluster_by_brain_pct_melt,aes(x = variable,y = value,fill = CellType)) +
  geom_bar(position = "stack",stat = "identity") +
  scale_fill_manual(values = cluster_cols) +
  labs(x = "Donor",
       y = "Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(plot = cluster_brain_bar,
       filename = here("plots","12_snRNA","Dim_Red","CellType_Final_By_BrainID_Ordered.pdf"),
       height = 12,
       width = 12)


#Violin plots by % mito, # genes, # UMIs, doubletScore
sce$CellType.Final.rev <- factor(x = sce$CellType.Final,
                                 levels = rev(c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                                                "DRD2_MSN_A","DRD2_MSN_B",
                                                "Inh_A","Inh_B","Inh_C","Inh_D","Inh_E","Inh_F",
                                                "Excitatory",
                                                "Astrocyte_A","Astrocyte_B","Ependymal",
                                                "Oligo","OPC",
                                                "Microglia",
                                                "Endothelial")))

# % mito
pct_mito_vln <- plotColData(object = sce,x = "subsets_Mito_percent",
                            y = "CellType.Final.rev",colour_by = "CellType.Final") +
  scale_color_manual(values = cluster_cols) +
  theme(legend.position = "none") +
  labs(y = "% Mitochondria",
       x = "Cell Type") +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)

ggsave(plot = pct_mito_vln,
       filename = here("plots","12_snRNA","Supplementary","Pct_Mito_by_CellType_Violin.pdf"))

# Number of genes
Number_Genes_vln <- plotColData(object = sce,x = "detected",
                                y = "CellType.Final.rev",colour_by = "CellType.Final") +
  scale_color_manual(values = cluster_cols) +
  scale_y_log10() +
  theme(legend.position = "none") +
  labs(y = "log10(Number of Genes)",
       x = "Cell Type") +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)

ggsave(plot = Number_Genes_vln,
       filename = here("plots","12_snRNA","Supplementary","Number_Genes_by_CellType_Violin.pdf"))

# Library Size
lib_size_vln <- plotColData(object = sce,x = "sum",
                                y = "CellType.Final.rev",colour_by = "CellType.Final") +
  scale_color_manual(values = cluster_cols) +
  scale_y_log10() +
  theme(legend.position = "none") +
  labs(y = "log10(Library Size)",
       x = "Cell Type") +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)

ggsave(plot = lib_size_vln,
       filename = here("plots","12_snRNA","Supplementary","Library_Size_by_CellType_Violin.pdf"))

# doubletScore
dub_Score_vln <- plotColData(object = sce,x = "doubletScore",
                            y = "CellType.Final.rev",colour_by = "CellType.Final") +
  scale_color_manual(values = cluster_cols) +
  theme(legend.position = "none") +
  labs(y = "doubletScore",
       x = "Cell Type") +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)

ggsave(plot = dub_Score_vln,
       filename = here("plots","12_snRNA","Supplementary","doubletScore_by_CellType_Violin.pdf"))


####Sort by celltype
#Calculate cluster percentage by sort
sort_by_cluster <- as.data.frame.matrix(table(sce$Sort,sce$CellType.Final))

#Calculate percentages 
sort_by_cluster_pct <- sweep(sort_by_cluster,MARGIN = 2,colSums(sort_by_cluster),"/") * 100
sort_by_cluster_pct$Sort <- rownames(sort_by_cluster_pct)

#Melt dataframe and plot
sort_by_cluster_pct_melt <- reshape2::melt(sort_by_cluster_pct)

#make some new sort colors
sort_cols <- Polychrome::createPalette(length(unique(sce$Sort)),
                                       c("#D81B60","#004D40"))
names(sort_cols) <- unique(sce$Sort)

sort_cluster_bar <- ggplot(data = sort_by_cluster_pct_melt,aes(x = variable,y = value,fill = Sort)) +
  geom_bar(position = "stack",stat = "identity") +
  scale_fill_manual(values = sort_cols) +
  labs(x = "CellType",
       y = "Percent") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = sort_cluster_bar,
       filename = here("plots","12_snRNA","Dim_Red","Sort_by_CellType_Final_Ordered.pdf"),
       height = 12,
       width = 12)

#Remake above plot but with sort on the x-axis and celltype proportion on y-axis
cluster_by_sort <- as.data.frame.matrix(table(sce$CellType.Final,sce$Sort))

#Calculate percentages 
cluster_by_sort_pct <- sweep(cluster_by_sort,MARGIN = 2,colSums(cluster_by_sort),"/") * 100
cluster_by_sort_pct$CellType <- rownames(cluster_by_sort_pct)

#Melt dataframe and plot
cluster_by_sort_pct_melt <- reshape2::melt(cluster_by_sort_pct)

#Make the bargraph.
sort_cluster_bar <- ggplot(data = cluster_by_sort_pct_melt,aes(x = variable,y = value,fill = CellType)) +
  geom_bar(position = "stack",stat = "identity") +
  scale_fill_manual(values = cluster_cols) +
  labs(x = "Sort",
       y = "Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(plot = sort_cluster_bar,
       filename = here("plots","12_snRNA","Dim_Red","CellType_By_Sort_Final_Ordered.pdf"),
       height = 12,
       width = 12)


####Make sex on x-axis and percentage of total cells 
#This is the percentage of each sex by cluster. 
#First need to make a column for Sex
sce$Sex <- ifelse(sce$Brain_ID %in% c("Br2720","Br8667","Br8492","Br8325"),
                  "Female",
                  "Male")

#Calculate cluster percentage by BrainID
sex_by_cluster <- as.data.frame.matrix(table(sce$CellType.Final,sce$Sex))

#Calculate percentages 
sex_by_cluster_pct <- sweep(sex_by_cluster,MARGIN = 2,STATS = colSums(sex_by_cluster),FUN = "/") * 100
sex_by_cluster_pct$CellType <- rownames(sex_by_cluster_pct)

#Melt dataframe and plot
sex_by_cluster_pct_melt <- reshape2::melt(sex_by_cluster_pct)
#Using CellType as id variables


sex_cluster_bar <- ggplot(data = sex_by_cluster_pct_melt,aes(x = variable,y = value,fill = CellType)) +
  geom_bar(position = "stack",stat = "identity") +
  scale_fill_manual(values = cluster_cols) +
  labs(x = "Sex",
       y = "Percent") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(plot = sex_cluster_bar,
       filename = here("plots","12_snRNA","Dim_Red","Sex_by_CellType_Final_Ordered.pdf"))



#Violin of general marker genes. 
gen_vln <- plotExpression(object = sce,
                          features = c("RBFOX3","GJA1","MOBP","PDGFRA","C3","EBF1"),
                          swap_rownames = "gene_name",
                          x = "CellType.Final", colour_by = "CellType.Final",
                          ncol = 1) +
  scale_color_manual(values = cluster_cols) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") 

ggsave(plot = gen_vln,
       filename = here("plots","12_snRNA","Expression","General_Violin_Figure1.pdf"),
       height = 7,width = 7)

ggsave(plot = gen_vln,
       filename = here("plots","12_snRNA","Expression","General_Violin_Figure1.png"),
       height = 7,width = 7)

### Read in the pairwise DEGs
load(here("processed-data","12_snRNA","markers_pairwise_CellType_Final.rda"),verbose = TRUE)
# Loading objects:
#   markers_pairwise

#Make CSVs for the marker genes. 
for(i in names(markers_pairwise)){
  print(i)
  x <- markers_pairwise[[i]]
  write.csv(x = x,
            file = here("processed-data","12_snRNA","DEG_CSVs","Pairwise",paste0(i,"_Pairwise.csv")))
}


### Read in the DEGs
load(here("processed-data","markers_1vAll_CellType_Final.rda"),verbose = TRUE)
# Loading objects:
#   markers_1vALL_df

#Make CSVs for the marker genes. 
for(i in unique(markers_1vALL_df$cellType.target)){
  print(i)
  x <- subset(markers_1vALL_df,subset=(cellType.target == i))
  write.csv(x = x,
            file = here("processed-data","12_snRNA","DEG_CSVs","1vALL",paste0(i,"_1vALL.csv")))
}

###########Complext heatmap of basic markers
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(sce$CellType.Final)

############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.  
markers_all <- c("SYT1","SNAP25","GAD1","PPP1R1B", #Broad neurons
                "DRD1","RXFP1", #D1/D1 islands
                "RELN","CNTNAP3B", #D1_A
                "TRHDE","CPNE4",#D1_B
                "RXFP2","SEMA3E", #D1_C
                "VWC2L","CLSTN2", #D1_D
                "DRD2","ADORA2A", #General D2
                "PENK","PTPRM",#D2_A
                "CLMP","GRIK3",#D2_B
                "IL1RAPL2","PDGFD", #Inh_A
                "VIP","CCK", #Inh_B
                "GLP1R","TAC3", #Inh_C
                "CHAT","SLC5A7", #Inh_D
                "NPY","SST", #Inh_E
                "KCNC2","ANK1", #Inh_F
                "SLC17A7","TBR1", #Excitatory
                "AQP4","GFAP", #Astrocyte
                "CAPS","FOXJ1", #Ependymal
                "MBP","MOBP", #Oligo
                "PDGFRA","VCAN", #OPC
                "C3","DOCK8", #Microglia
                "DCN","CLDN5") #Endothelial  

#marker labels
marker_labels <- c(rep("Neuron",4),
                   rep("D1_MSN",10),
                   rep("D2_MSN",6),
                   rep("Inhibitory",12),
                   rep("Excitatory",2),
                   rep("Astrocyte",2),
                   rep("Ependymal",2),
                   rep("Oligo",2),
                   rep("OPC",2),
                   rep("Microglia",2),
                   rep("Endothelial",2))

marker_labels <- factor(x = marker_labels,
                        levels =  unique(marker_labels))

colors_markers <- list(marker = c(Neuron = "black",
                                  D1_MSN = "#332288",
                                  D2_MSN = "#D81B60",
                                  Inhibitory = "#44AA99",
                                  Excitatory = as.character(cluster_cols["Excitatory"]),
                                  Astrocyte = "#DDCC77",
                                  Ependymal = as.character(cluster_cols["Ependymal"]),
                                  Oligo = as.character(cluster_cols["Oligo"]),
                                  OPC = as.character(cluster_cols["OPC"]),
                                  Microglia = as.character(cluster_cols["Microglia"]),
                                  Endothelial = as.character(cluster_cols["Endothelial"])))

col_ha <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                           show_annotation_name = FALSE,
                                           show_legend = TRUE,
                                           col = colors_markers)

###########set up rows for heatmap. 
# cluster labels
cluster_pops <- list(D1_MSN = c("DRD1_MSN_A","DRD1_MSN_B",
                                "DRD1_MSN_C","DRD1_MSN_D"),
                     D2_MSN = c("DRD2_MSN_A","DRD2_MSN_B"),
                     Inhibitory = c("Inh_A","Inh_B","Inh_C",
                                    "Inh_D","Inh_E","Inh_F"),
                     Excitatory = "Excitatory",
                     Astrocyte = c("Astrocyte_A","Astrocyte_B"),
                     Ependymal = "Ependymal",
                     Oligo = "Oligo",
                     OPC = "OPC",
                     Microglia = "Microglia",
                     Endothelial = "Endothelial")

# cluster labels order
# # cluster labels order
cluster_pops_order <- unname(unlist(cluster_pops))

# swap values and names of list
cluster_pops_rev <- rep(names(cluster_pops),
                        times = sapply(cluster_pops, length))
names(cluster_pops_rev) <- unname(unlist(cluster_pops))
#cluster_pops_rev <- cluster_pops_rev[as.character(sort(cluster_pops_order))]
cluster_pops_rev <- factor(cluster_pops_rev, levels = names(cluster_pops))

# second set of cluster labels
neuron_pops <- ifelse(cluster_pops_rev %in% c("Inhibitory","D1_MSN",
                                              "D2_MSN","Excitatory"),
                      "Neuron",
                      "Non-neuron")

neuron_pops <- factor(x = neuron_pops,levels = c("Neuron","Non-neuron"))

colors_neurons <- list(class = c(Neuron = "black",
                                 `Non-neuron` = "gray65"))

#n <- table(sce$CellType.Final)

#row annotation dataframe. 
# row annotation
pop_markers <- list(population = c(D1_MSN = "#332288",
                                   D2_MSN = "#D81B60",
                                   Inhibitory = "#44AA99",
                                   Excitatory = as.character(cluster_cols["Excitatory"]),
                                   Astrocyte = "#DDCC77",
                                   Ependymal = as.character(cluster_cols["Ependymal"]),
                                   Oligo = as.character(cluster_cols["Oligo"]),
                                   OPC = as.character(cluster_cols["OPC"]),
                                   Microglia = as.character(cluster_cols["Microglia"]),
                                   Endothelial = as.character(cluster_cols["Endothelial"])))

row_ha <- rowAnnotation(class = neuron_pops,
                        population = cluster_pops_rev,
                        show_annotation_name = FALSE,
                        show_legend = FALSE,
                        col = c(pop_markers,colors_neurons))


dat <- assay(sce,"logcounts")
rownames(dat) <- rowData(sce)$gene_name
dim(dat)
#[1]  36601 103339

dat <- dat[markers_all,]
dim(dat)
#[1]     46 103339

dat <- as.matrix(dat)

hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))),
                center = TRUE,
                scale = TRUE)

hm_mat <- hm_mat[names(cluster_pops_rev),]


min(hm_mat)
#[1] -1.478196

max(hm_mat)
#[1] 4.248515

col_fun <- circlize::colorRamp2(c(-2,0,4),c("blue","white","red"))

hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "centered,scaled",
                              column_title = "General cell class marker \ngene expression across clusters",
                              column_title_gp = gpar(fontface = "bold"),
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              bottom_annotation = col_ha,
                              right_annotation = row_ha,
                              column_split = marker_labels,
                              row_split = cluster_pops_rev,
                              row_title = NULL,
                              rect_gp = gpar(col = "gray50", lwd = 0.5),
                              col = col_fun)


pdf(file = here("plots","12_snRNA","Expression","CellType_Final_Heatmap_BroadMarkers.pdf"),
    width = 16,
    height = 12)
draw(hm)
dev.off()


###Feature plots for SNAP25, GAD1, DRD1, and DRD2
#snap25
snap25_featureplot <- plotReducedDim(sce,
                                     dimred = "tSNE_HARMONY",
                                     colour_by = "SNAP25",
                                     swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","red"))

ggsave(plot = D1_featureplot,
       filename = here("plots","12_snRNA","Expression","DRD1_Featureplot_Figure1.pdf"),
       height = 8,width = 8)

#DRD1
D1_featureplot <- plotReducedDim(sce,
                                 dimred = "tSNE_HARMONY",
                                 colour_by = "DRD1",
                                 swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","red"))

ggsave(plot = D1_featureplot,
       filename = here("plots","12_snRNA","Expression","DRD1_Featureplot_Figure1.pdf"),
       height = 8,width = 8)

#DRD2
D2_featureplot <- plotReducedDim(sce,
                                 dimred = "tSNE_HARMONY",
                                 colour_by = "DRD2",
                                 swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","red"))

ggsave(plot = D2_featureplot,
       filename = here("plots","12_snRNA","Expression","DRD2_Featureplot_Figure1.pdf"),
       height = 8,width = 8)

