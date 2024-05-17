#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(dplyr)
library(scran)
library(here)

##############
#Load sce
print("Loading object")
Sys.time()
sce <- readRDS(here("processed-data","12_snRNA","sce_clustered_log.Rds"))
Sys.time()

dim(sce)

identical(rownames(colData(sce)),colnames(sce))

sce

########
#Preannotation plots
#Calculate cluster percentage by BrainID
brain_by_cluster <- as.data.frame.matrix(table(sce$Brain_ID,sce$k_20_walktrap))

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
  labs(x = "k=20 walktrap cluster",
       y = "Percent") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = brain_cluster_bar,
       filename = here("plots","12_snRNA","Dim_Red","BrainID_by_cluster_k_20_walktrap_preannotation.png"),
       height = 12,
       width = 12)

#Number of cells barplot
Num_Nucs <- table(sce$k_20_walktrap) %>% as.data.frame() %>% 
  ggplot(aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "k=20 walktrap cluster",
       y = "Number of Nuclei") 

ggsave(plot = Num_Nucs,
       filename = here("plots","12_snRNA","k_20_walktrap_NumberofNuclei_Bargraph.png"),
       height = 12,
       width = 12)

#Number of cells barplot split by sort
Num_Nucs_sort <- table(sce$k_20_walktrap,sce$Sort) %>% as.data.frame() %>% 
  ggplot(aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "k=20 walktrap cluster",
       y = "Number of Nuclei") 
ggsave(plot = Num_Nucs_sort,
       filename = here("plots","12_snRNA","k_20_walktrap_NumberofNucleibySort_Bargraph.png"),
       height = 12,
       width = 12)

#Number of genes detected by cluster
detected_vln <- plotColData(object = sce,
                            y = "detected",
                            x = "k_20_walktrap",
                            colour_by = "k_20_walktrap") +
  ggtitle("Number of Detected Genes") +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) 
ggsave(plot = detected_vln,
       filename = here("plots","12_snRNA","k_20_walktrap_detected.png"),
       height = 12,
       width = 12)

#percent mito
percent_mito_vln <- plotColData(object = sce,
                                y = "subsets_Mito_percent",
                                x = "k_20_walktrap",
                                colour_by = "k_20_walktrap") +
  ggtitle("Percent mito") +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) 
ggsave(plot = percent_mito_vln,
       filename = here("plots","12_snRNA","k_20_walktrap_PercentMito.png"),
       height = 12,
       width = 12)

######################################
#Cluster Annotations
#read in the annotation frame
annotation_df <- read.csv(file = here("processed-data","12_snRNA",
                                      "Annotation_dataframe.csv"),
                          header = TRUE)


sce$CellType.Final <- annotation_df$Annotation[match(sce$k_20_walktrap,
                                                     annotation_df$Cluster)]


#Annotate the tSNE
cluster_cols <- Polychrome::createPalette(length(unique(sce$CellType.Final)),
                                          c("#D81B60", "#1E88E5","#FFC107","#009E73"))
names(cluster_cols) <- unique(sce$CellType.Final)

save(cluster_cols,file = here("processed-data","12_snRNA","CellType_Final_Cluster_Cols.rda"))

#tSNE with annotations
annotated_tSNE <- plotReducedDim(object = sce,
                                 dimred = "tSNE_HARMONY",
                                 colour_by = "CellType.Final",
                                 text_by = "CellType.Final") +
  scale_color_manual(values = cluster_cols) +
  theme(legend.position = "none")
ggsave(filename = here("plots","12_snRNA","tSNE_CellType_Final.png"),
       plot = annotated_tSNE,height = 12, width = 12)

#umap with annotations
annotated_umap <- plotReducedDim(object = sce,
                                 dimred = "umap_HARMONY",
                                 colour_by = "CellType.Final",
                                 text_by = "CellType.Final") +
  scale_color_manual(values = cluster_cols) +
  theme(legend.position = "none")
ggsave(filename = here("plots","12_snRNA","umap_CellType_Final.png"),
       plot = annotated_umap,height = 12, width = 12)

##########################3
#Post annotation plots
#Calculate cluster percentage by BrainID
brain_by_cluster <- as.data.frame.matrix(table(sce$Brain_ID,sce$CellType.Final))

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
       filename = here("plots","12_snRNA","Dim_Red","BrainID_by_CellType_Final.png"),
       height = 12,
       width = 12)

#Number of cells barplot
Num_Nucs <- table(sce$CellType.Final) %>% as.data.frame() %>% 
  ggplot(aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_cols) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell Type",
       y = "Number of Nuclei") 

ggsave(plot = Num_Nucs,
       filename = here("plots","12_snRNA","CellType.Final_NumberofNuclei_Bargraph.png"),
       height = 12,
       width = 12)

#Number of cells barplot split by sort
Num_Nucs_sort <- table(sce$CellType.Final,sce$Sort) %>% as.data.frame() %>% 
  ggplot(aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell Type",
       y = "Number of Nuclei") 
ggsave(plot = Num_Nucs_sort,
       filename = here("plots","12_snRNA","CellType.Final_NumberofNucleibySort_Bargraph.png"),
       height = 12,
       width = 12)

#Number of genes detected by cluster
detected_vln <- plotColData(object = sce,
                            y = "detected",
                            x = "CellType.Final",
                            colour_by = "CellType.Final") +
  ggtitle("Number of Detected Genes") +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) +
  scale_color_manual(values = cluster_cols) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) 
ggsave(plot = detected_vln,
       filename = here("plots","12_snRNA","CellType_Final_detected.png"),
       height = 12,
       width = 12)

#percent mito
percent_mito_vln <- plotColData(object = sce,
                                y = "subsets_Mito_percent",
                                x = "CellType.Final",
                                colour_by = "CellType.Final") +
  ggtitle("Percent mito") +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) + 
  scale_color_manual(values = cluster_cols) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) 
ggsave(plot = percent_mito_vln,
       filename = here("plots","12_snRNA","CellType_Final_PercentMito.png"),
       height = 12,
       width = 12)

#####
#save the object with cell type information
print("Saving Object with residuals included")
Sys.time()
saveRDS(sce,file = here("processed-data","12_snRNA","sce_CellType.Rds"))
Sys.time()
print("Object saved")

#Now remove the binomial_pearson_residuals and save a separate object
assay(sce,"binomial_pearson_residuals") <- NULL
print("Saving Object without residuals")
Sys.time()
saveRDS(sce,file = here("processed-data","12_snRNA","sce_CellType_noresiduals.Rds"))
Sys.time()
print("Object saved")

#Repro information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
