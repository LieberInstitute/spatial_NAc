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

#calculate purity 
pure.sce <- neighborPurity(reducedDim(sce, "HARMONY"), sce$k_10_louvain_1)

#configure dataframe
pure.data <- as.data.frame(pure.sce)
pure.data$maximum <- factor(pure.data$maximum)
pure.data$cluster <- sce$k_10_louvain_1

purity_plot <- ggplot(pure.data, aes(x=cluster, y=purity, colour=maximum)) +
                ggbeeswarm::geom_quasirandom(method="smiley")

ggsave(plot = purity_plot,file = here("plots","12_snRNA","Dim_Red","Purity_beeswarm.png"),height = 12, width = 12)

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
