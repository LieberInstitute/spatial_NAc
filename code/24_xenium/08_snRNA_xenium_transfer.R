#Perform label transfer between the snRNA-seq data and the seurat objects for each donor
library(Seurat)
library(here)

#Raise limit
options(future.globals.maxSize = 16 * 1024^3) #16 GB

#load snRNA-seq seurat object
snRNA <- readRDS(here("processed-data", "12_snRNA","seurat_snrna_cellType_noResiduals.Rds"))

snRNA

#load the annotated seurat object containing banksy clusters first. 
xen <- readRDS(here("processed-data", "xenium", 
                    "SPEs", "seurat_xenium.Rds"))

#Perform label transfer
#Will use top 30 PCs for initial trial of Seurat label transfer
#Identify the transfer anchors
message(paste0("Finding transfer anchors - ", Sys.time()))
anchors <- FindTransferAnchors(reference = snRNA, 
                               query = xen, 
                               dims = 1:30,
                               reference.reduction = "seurat_pca",
                               reduction = "rpca")

#Perform the label transfer
message(paste0("Transferring Data - ", Sys.time()))
predictions <- TransferData(anchorset = anchors, 
                            refdata = snRNA$CellType.Final, 
                            dims = 1:30)

#Add predictions to the object
message(paste0("Adding metadata - ", Sys.time()))
xen <- AddMetaData(xen, 
                   metadata = predictions)

#Save the predictions  and updated seurat objectas an RDS file. 
message(paste0("Saving objects - ", Sys.time()))
saveRDS(object = predictions,file = here("processed-data","12_snRNA","xenium_snRNA_predictions.Rds"))
saveRDS(object = xen,file = here("processed-data", "12_snRNA","xenium_snRNA_predictions_object.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
