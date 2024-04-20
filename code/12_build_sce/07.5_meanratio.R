#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac
library(SingleCellExperiment)
library(DeconvoBuddies)
library(HDF5Array)
library(here)

#Load sce
print("Loading SCE")
sce <- loadHDF5SummarizedExperiment(dir = here("processed-data","12_snRNA",
                                               "sce_clustered"))
dim(sce)

stopifnot(identical(rownames(colData(sce)),colnames(sce)))

sce
######################

#Run mean ratio
gmr_pt25 <- get_mean_ratio2(sce,
                            cellType_col = "k_10_louvain_pt25",
                            assay_name = "logcounts",
                            add_symbol = FALSE)

save(gmr_pt25,
     file = here("processed-data","12_snRNA","broad_clustering_meanratio_preannotation.rds"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
