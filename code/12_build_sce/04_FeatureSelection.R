#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

library(SingleCellExperiment)
library(sessioninfo)
library(HDF5Array)
library(ggplot2)
library(scater)
library(scran)
library(scry)
library(here)

#Load the QCed sce object
sce <- loadHDF5SummarizedExperiment(dir = here("processed-data","12_snRNA",
                                               "sce_clean"))

dim(sce)

#Deviance feature selection
sce <- devianceFeatureSelection(sce,
                                assay = "counts",
                                fam = "binomial",
                                sorted = FALSE,
                                batch = as.factor(sce$Sample))

#Plot binomial deviance
pdf(here("plots","12_snRNA","featureSelxn_binomialDeviance-byGene_numericcutoffs.pdf"))
plot(sort(rowData(sce)$binomial_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylab = "binomial deviance"
)
abline(v = 2000,lty = 2, col = "red")
dev.off()


#Calculate null residuals
sce <- nullResiduals(sce,
                     assay = "counts", 
                     fam   = "binomial", 
                     type  = "deviance")

#Save object
saveHDF5SummarizedExperiment(sce,
                             here("processed-data","12_snRNA",
                                  "sce_featureselection"),
                             replace = TRUE)


print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

