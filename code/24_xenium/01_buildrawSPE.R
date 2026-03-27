#Goal: Build raw SPE
# cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
# module load r_nac
#code modified from https://github.com/LieberInstitute/spatialDLPFC_SCZ_XENIUM/blob/devel/code/analysis/01_build_spe/01_build_spe.R

library(SingleCellExperiment)
library(SpatialExperiment)
library(sessioninfo)
library(tidyverse)
library(escheR)
library(here)

#Load sample_info dataframe
sample_info <- read.csv(here("processed-data","xenium","Spatial_NAc_Xenium_Sample_Info.csv"))

#Add the full sample path to the sample_info dataframe
sample_info$full_data_path <- file.path(paste(here("raw-data","xenium"),sample_info$Output_Directory,sep = "/"))

all_spes <- vector(mode = "list")
names(all_spes) <- sample_info$Sample_ID
for(i in 1:nrow(sample_info)){
  print(i)
  #Pull path to data and sample_ID
  sample_path <- sample_info[i,"full_data_path"]
  Sample <- sample_info[i,"Sample"]
  
  print(sprintf("----------%s----------", Sample))
  
  sample_info_use <- sample_info[i,]
  
  #Paths to counts matrix and cell csv file
  counts_path <- here(sample_path,"cell_feature_matrix.h5")
  cell_info_path <- here(sample_path,"cells.csv.gz")
  
  #Creae sce object
  sce <- DropletUtils::read10xCounts(samples = counts_path, sample.names = Sample)
  counts(sce) <- methods::as(DelayedArray::realize(counts(sce)), "dgCMatrix") # Convert to delayed array
  
  cell_info <- vroom::vroom(cell_info_path) #reads in cell csv file much faster than standard functions
  
  #Add cell info to sce
  colData(sce) <- cbind(colData(sce),cell_info)
  spe <- toSpatialExperiment(sce,spatialCoordsNames = c("x_centroid","y_centroid"))
  rownames(spe) <- rowData(spe)$Symbol #No need to uniquify because Ensembl speicifed during panel generation
  
  #Need to assign donor-specific column names because we will merge everything after loop is finished. 
  colnames(spe) <- paste(Sample,rownames(cell_info),sep = "_")
  
  #Add in relevant metadata. 
  #First technical variables. 
  spe$Xenium_Run_ID <- sample_info_use$Xenium_Run_ID
  spe$Slide_ID <- sample_info_use$Slide_ID
  
  #Now biological variables
  #spe$sample_id <- Sample
  spe$Age <- sample_info_use$AgeDeath
  spe$Sex <- sample_info_use$Sex
  spe$PrimaryDx <- sample_info_use$PrimaryDx
  
  #Plot a pre-rotation image. 
  x <- escheR::make_escheR(spe) |>
    escheR::add_fill("total_counts") +
    scale_fill_gradient(low = "white", high = "black")
  ggsave(plot = x,
         filename = here("plots","xenium","01_buildSPE",
                         paste0(Sample,"_prerotation.png")),
         height = 16, width = 18)
  
  #Pull spatial coordinates
  coords <- spatialCoords(spe)
  
  #Pull the center of x and y coordinates
  center_x <- mean(range(coords[,1]))  
  center_y <- mean(range(coords[,2]))
  
  #Generate final coordinates
  x_final <- center_x - (coords[,1] - center_x) #Identify distance of each x from center and subtract from center to mirror horizontally. 
  y_final <- center_y - (coords[,2] - center_y) #Identify distance of each y from center and subtract from center to mirror vertically. 
  
  #Convert the coordinates
  spatialCoords(spe) <- cbind(x_final,y_final)
  
  #Plot total_counts on top of the tissue to visualize rotation. 
  x <- escheR::make_escheR(spe) |>
    escheR::add_fill("total_counts") +
    scale_fill_gradient(low = "white", high = "black")
  ggsave(plot = x,
         filename = here("plots","xenium","01_buildSPE",
                         paste0(Sample,"_rotated.png")),
         height = 16, width = 18)
  
  
  all_spes[[Sample]] <- spe
  
  #clear memory
  rm(sce)
  rm(cell_info)
  gc() #garbage collection
}

#Combine all of the spes 
message(paste0("Combining the spes - ",Sys.time()))
all_spes <- do.call(cbind,all_spes)

message(paste0("Saving SPE - ",Sys.time()))
saveRDS(all_spes,here("processed-data","xenium","SPEs","spe_raw.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
