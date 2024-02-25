#Concatenate samples and build basic sce
#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc
#module load r_nac

library(SingleCellExperiment)
library(sessioninfo)
library(DropletUtils)
library(rtracklayer)
library(here)

#Read in a dataframe consisting of identifying information for the samples
sample_data <- read.delim(here("processed-data",
                               "12_snRNA",
                               "Spatial_NAc_snRNA_Seq_sample_info.csv"),
                          header = TRUE,sep = ",")

#Add raw data paths to the sample dataframe
sample_data$Raw_data_path <- NA
for(i in 1:20){
    sample_data[i,"Raw_data_path"] <-paste0("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/11_cellranger/",
                                            i,
                                            "c_NAc_SVB/outs/raw_feature_bc_matrix")
}

#Check that all files exist. Also, stop if they do not. 
all(file.exists(sample_data$Raw_data_path))

stopifnot(all(file.exists(sample_data$Raw_data_path)))

## Build basic SCE
message("Read 10x data and create sce - ", Sys.time())
sce <- read10xCounts(samples = sample_data$Raw_data_path,
                     sample.names = sample_data$Sample,
                     type = "sparse",
                     col.names = TRUE)
message("RDone - ", Sys.time())

######colData
#merging removes the rownames that are unique to each sample/cell. 
#create a key that can be used to reorder the column data to the original order that it was read in
sce$key <- paste(sce$Barcode,sce$Sample,sep = "_")

#Add relevant information to the columnData
#Key raw data path information intact. Additional sanity check for later. 
new_column_data<- merge(x  = colData(sce),
                        y  = sample_data,
                        by = "Sample")

#reorder the colData so that it is in the original order. 
new_column_data <- new_column_data[match(sce$key, new_column_data$key), ]

#Check that the key column within the new column data dataframe and sce are in the same order
#Stop if not
identical(sce$key,new_column_data$key)

stopifnot(identical(sce$key, new_column_data$key))

#Add the rownames back.
#rownames for the sce object need to be sample#_Barcode 
rownames(new_column_data) <- paste(as.character(lapply(strsplit(new_column_data$Sample,split = "c"),"[",1)),new_column_data$Barcode,sep = "_")

#Update the column data. 
colData(sce) <- new_column_data

#rownames of the coldata needs to be in same order as column names of the count matrix
identical(rownames(colData(sce)),colnames(sce))

stopifnot(identical(rownames(colData(sce)),colnames(sce)))

####rowData
gtf <- rtracklayer::import("/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf")
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

#match the genes
match_genes <- match(rownames(sce),gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

#Keep only specific columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")]

#Add gene info
rowRanges(sce) <- gtf[match_genes]

#Check out object so far
sce

#Empty droplets have not been removed. Will be the next step of the analysis. 

#Save object. 
save(sce,
     file = here("processed-data","12_snRNA","sce_raw.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
