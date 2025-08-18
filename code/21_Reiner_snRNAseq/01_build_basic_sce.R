library(SingleCellExperiment)
library(Seurat)
library(sessioninfo)
library(DropletUtils)
library(rtracklayer)
library(here)

# Specify input directories
rawDir <- "/dcs05/lieber/marmaypag/temp_data/01-nac/processed-data/01-cellranger"
samples_list <- list.dirs(rawDir,recursive = FALSE)
samples_list <- gsub(paste0(rawDir, "/"),"", samples_list)

samples_df <- data.frame(samples_list)
colnames(samples_df) <- "sample"
samples_df$Raw_data_path <- NA
for(i in 1:dim(samples_df)[1]){
    samples_df[i,"Raw_data_path"] <-paste0(rawDir, "/", samples_df$sample[i],
                                            "/outs/raw_feature_bc_matrix")
}

#Check that all files exist. Also, stop if they do not. 
all(file.exists(samples_df$Raw_data_path))
stopifnot(all(file.exists(samples_df$Raw_data_path)))

## Build basic SCE
message("Read 10x data and create sce - ", Sys.time())
sce <- read10xCounts(samples = samples_df$Raw_data_path,
                     sample.names = samples_df$sample,
                     type = "sparse",
                     col.names = TRUE)
message("RDone - ", Sys.time())

# Read in gene data and add to this object
####rowData
gtf <- rtracklayer::import("/dcs05/lieber/marmaypag/temp_data/01-nac/refs/refdata-gex-mRatBN7-2-2024-A/genes/genes.gtf")
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

#match the genes
match_genes <- match(rownames(sce),gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

#Keep only specific columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_biotype")]

#Add gene info
rowRanges(sce) <- gtf[match_genes]
#Save object. 
saveRDS(sce,
        file = here("processed-data","21_Reiner_snRNAseq","sce_raw.Rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
