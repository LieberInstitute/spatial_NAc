#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac
library(SingleCellExperiment)
library(DeconvoBuddies)
library(here)

#Load sce
print("Loading SCE")
Sys.time()
sce <- readRDS(here("processed-data","12_snRNA","sce_CellType_noresiduals.Rds"))
dim(sce)

stopifnot(identical(rownames(colData(sce)),colnames(sce)))

sce
######################

#Run mean ratio
gmr_CT <- get_mean_ratio2(sce,
                          cellType_col = "CellType.Final",
                          assay_name = "logcounts",
                          add_symbol = FALSE)

save(gmr_CT,
     file = here("processed-data","12_snRNA","CellType_Final_meanratio.rda"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
