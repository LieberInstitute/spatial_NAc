# cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
# module load r_nac

library(SingleCellExperiment)
library(pheatmap)
library(bluster)
library(scater)
library(scran)
library(here)

#Load in the object
print("Loading Object")
Sys.time()
sce <- readRDS(here("processed-data","12_snRNA","sce_clustered_log.Rds"))
print("Object loaded")
Sys.time()

sce

#Build the graph
Sys.time()
set.seed(10)
snn_k_10 <- buildSNNGraph(sce, k = 10, use.dimred = "HARMONY")
Sys.time()

# "Rather, we use the pairwiseModularity() function from bluster with as.ratio=TRUE, 
# which returns the ratio of the observed to expected sum of weights between each pair of clusters. 
# We use the ratio instead of the difference as the former is less dependent on the number of cells 
# in each cluster" - OSCA advanced, Chatper 5.2.5
k_10_modularity <- bluster::pairwiseModularity(graph = snn_k_10,
                                               clusters = sce$k_10_louvain_1,
                                               as.ratio = TRUE)


save(k_10_modularity,file = here("processed-data","12_snRNA","k_10_modularity_scores.rda"))

##Generate the heatmap
pdf(file = here("plots","12_snRNA","Dim_Red","k_10_pairwise_modularity_final_celltypes.pdf"),height = 15, width =15)
pheatmap(log2(k_10_modularity+1), 
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         display_numbers=TRUE, 
         number_format="%.2f", 
         fontsize_number=6.5,
         main = "Modularity ratio for 27 graph-based clusters in human NAc (n=20)",
         color=colorRampPalette(c("white","orange","red"))(100))
dev.off()


print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
