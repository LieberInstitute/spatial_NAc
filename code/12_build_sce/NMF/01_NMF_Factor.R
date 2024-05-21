# cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/12_snRNA/code/NMF
# module load r_nac

#Load libraries
library(SingleCellExperiment)
library(here)
library(scuttle)
library(Matrix)
library(sessioninfo)
library("RcppML",lib.loc = "/users/rphillip/R/4.3.x") #Need development version for singlet. 
library("singlet",lib.loc = "/users/rphillip/R/4.3.x")

#load the object
print("Loading Object")
Sys.time()
sce <- readRDS(here("processed-data","12_snRNA","sce_CellType_noresiduals.Rds"))

sce

#Run cross validation
print("Running NMF")
Sys.time()
options(RcppML.threads=4)
x <- RcppML::nmf(assay(sce,"logcounts"),
                 k=66,
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 1234,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE)



#Save the cross validation results
saveRDS(x, file = here("processed-data","12_snRNA","NMF","NMF_Results.Rds"))

#sesion info
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
