# cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
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
print("Running cross validation")
Sys.time()
cvnmf <- cross_validate_nmf(
    logcounts(sce),
    ranks=c(10,50,100,125,150),
    n_replicates = 3,
    tol = 1e-03,
    maxit = 100,
    verbose = 3,
    L1 = 0.1,
    L2 = 0,
    threads = 0,
    test_density = 0.2
)

#Save the cross validation results
saveRDS(cvnmf, file = here("processed-data","12_snRNA","NMF","nmf_cv_results.Rds"))

#plot the results
pdf(here("plots","12_snRNA","NMF","cross_validation_results.pdf"))
plot(cvnmf)
dev.off()

#sesion info
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
