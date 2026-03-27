# cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
# module load r_nac
library(SingleCellExperiment)
library(ggplot2)
library(scater)
library(scran)
library(lisi)
library(here)


#Load the object 
#load the sce object
sce <- readRDS(here("processed-data","12_snRNA","sce_CellType_noresiduals.Rds"))

sce

#Remove neuron ambig
sce <- sce[,sce$CellType.Final != "Neuron_Ambig"]

embeddings_harmony <- reducedDim(sce,"HARMONY")[,1:2]
embeddings_tSNE <- reducedDim(sce,"tSNE_HARMONY")[,1:2]

harmony_lisi_brainID <- compute_lisi(embeddings_harmony, colData(sce), "Brain_ID",perplexity = 60)
harmony_lisi_Sample <- compute_lisi(embeddings_harmony, colData(sce), "Sample",perplexity = 60)

#Add lisi scores to the colData
sce$lisi_brainID <- harmony_lisi_brainID
sce$lisi_Sample <- harmony_lisi_Sample

#Make a ridge plot for the lisi scores
library(ggridges)
x <- as.data.frame(colData(sce))

#Load cluster colors 
load(here("processed-data","12_snRNA","070924_21colors_celltypeFinal.rda"),verbose = TRUE)

#Remove neuron_ambig
cluster_cols <- cluster_cols[-14]

#LISI ridge plots
#By BrainID
LISI_BrainID <- ggplot(data = x,aes(x = Brain_ID.1,y = CellType.Final,fill = CellType.Final)) +
  geom_density_ridges() +
  scale_fill_manual(values = cluster_cols) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none") +
  labs(x = "LISI",
       y = "Cell Type",
       title = "LISI by Brain ID\n10 Donors") 
ggsave(plot = LISI_BrainID,
       filename = here("plots","12_snRNA","Dim_Red","LISI_CellType_BrainID.png"),
       height = 8,width = 8)

#By Sample
LISI_Sample <- ggplot(data = x,aes(x = Sample.1,y = CellType.Final,fill = CellType.Final)) +
  geom_density_ridges() +
  scale_fill_manual(values = cluster_cols) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none") +
  labs(x = "LISI",
       y = "Cell Type",
       title = "LISI by Sample\n20 Samples") 
ggsave(plot = LISI_Sample,
       filename = here("plots","12_snRNA","Dim_Red","LISI_CellType_Sample.png"),
       height = 8,width = 8)


#Facet wrap the umap
umap_wrap <- as.data.frame(reducedDim(sce,"tSNE_HARMONY"))
umap_wrap$CellType <- sce$CellType.Final
umap_wrap$Sample <- sce$Sample
umap_wrap$Sort <- sce$Sort

sort_colors <- c("orange","forestgreen")
names(sort_colors) <- c("PI_NeuN","PI")

umap_wrap <- ggplot(data = umap_wrap,aes(x = TSNE1,y = TSNE2,color = Sort)) +
  geom_point(size = 0.75, alpha = 0.25) +
  scale_color_manual(values = sort_colors) +
  facet_wrap(~Sample) +
  theme_bw()
ggsave(plot = umap_wrap,
       filename = here("plots","12_snRNA","Dim_Red","tSNE_HARMONY_facetWrap.pdf"),
       height = 16,width = 16)

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
