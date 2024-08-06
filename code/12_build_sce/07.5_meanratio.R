#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac
library(SingleCellExperiment)
library(DeconvoBuddies)
library(here)

#Load sce
print("Loading SCE")
Sys.time()
sce <- readRDS(here("processed-data","12_snRNA","sce_clustered_log.Rds"))
dim(sce)

stopifnot(identical(rownames(colData(sce)),colnames(sce)))

sce
######################

#Run mean ratio
gmr_1 <- get_mean_ratio2(sce,
                         cellType_col = "k_10_louvain_1",
                         assay_name = "logcounts",
                         add_symbol = FALSE)

save(gmr_1,
     file = here("processed-data","12_snRNA","broad_clustering_meanratio_preannotation_k10_louvain_1.rda"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
