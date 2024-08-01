#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(scry)
library(here)

#Load the QCed sce object
sce <- readRDS(here("processed-data","12_snRNA","sce_clean.Rds"))

dim(sce)

class(sce$Sample)

sce

# feature selection
set.seed(717)
sce <- devianceFeatureSelection(sce,
                                assay = "counts",
                                fam = "binomial",
                                sorted = FALSE,
                                batch = sce$Sample)

#Plot binomial deviance
pdf(here("plots","12_snRNA","featureSelxn_binomialDeviance-byGene_numericcutoffs.pdf"))
plot(sort(rowData(sce)$binomial_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylab = "binomial deviance"
)
abline(v = 2000,lty = 2, col = "red")
dev.off()

#Calculate null residuals
#Pearson residuals. 
sce <- nullResiduals(sce,
                     assay = "counts", 
                     fam   = "binomial", 
                     type  = "pearson")

size_sce <- object.size(sce)
print(size_sce, units = "auto", standard = "IEC") 

#Save object
saveRDS(sce,
        file = here("processed-data","12_snRNA","sce_featureselection.Rds"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
