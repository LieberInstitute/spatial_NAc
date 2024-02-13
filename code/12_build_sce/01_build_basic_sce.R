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

## Build basic SCE
message("Read 10x data and create sce - ", Sys.time())
#Read 10x data and create sce - 2024-02-13 12:10:41.408777
sce <- read10xCounts(samples = sample_data$Raw_data_path,
                     sample.names = sample_data$Sample,
                     type = "sparse",
                     col.names = TRUE)
message("RDone - ", Sys.time())
#RDone - 2024-02-13 12:22:34.73507

######colData
#merging removes the rownames that are unique to each sample/cell. 
#create a column that will remain after merging. 
colData(sce)$unique_rowname <- rownames(colData(sce))

#Then merge and restore the rownames. 
new_column_data<- merge(x  = colData(sce),
                        y  = sample_data[,-which(colnames(sample_data) == "Raw_data_path")],
                        by = "Sample")
rownames(new_column_data) <- new_column_data$unique_rowname

#Update the column data. 
colData(sce) <- new_column_data


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
# class: SingleCellExperiment 
# dim: 36601 28269956 
# metadata(1): Samples
# assays(1): counts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(6): source type ... gene_name gene_type
# colnames(28269956): 10_AAACCCAAGAAACACT-1 10_AAACCCAAGAAACCCA-1 ...
# 9_TTTGTTGTCTTGGAAC-1 9_TTTGTTGTCTTTGGAG-1
# colData names(10): Sample Barcode ... Chromium_kit Library_Kit
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

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
# [1] "Reproducibility information:"
# [1] "2024-02-13 13:22:39 EST"
#     user   system  elapsed 
# 1848.925   59.797 4351.390 
# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-02-13
# pandoc   3.1.3 @ /jhpce/shared/libd/core/r_nac/1.0/nac_env/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.2)
# beachmat               2.18.0    2023-10-24 [1] Bioconductor
# Biobase              * 2.62.0    2023-10-24 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-01 [1] Bioconductor
# BiocIO                 1.12.0    2023-10-24 [1] Bioconductor
# BiocParallel           1.36.0    2023-10-24 [1] Bioconductor
# Biostrings             2.70.1    2023-10-25 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.2)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.2)
# codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.0)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# DelayedArray           0.28.0    2023-10-24 [1] Bioconductor
# DelayedMatrixStats     1.24.0    2023-10-24 [1] Bioconductor
# dqrng                  0.3.2     2023-11-29 [1] CRAN (R 4.3.2)
# DropletUtils         * 1.22.0    2023-10-24 [1] Bioconductor
# edgeR                  4.0.3     2023-12-10 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDb         * 1.38.1    2023-11-08 [1] Bioconductor
# GenomeInfoDbData       1.2.11    2023-12-12 [1] Bioconductor
# GenomicAlignments      1.38.0    2023-10-24 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-29 [1] Bioconductor
# HDF5Array              1.30.0    2023-10-24 [1] Bioconductor
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.2)
# IRanges              * 2.36.0    2023-10-24 [1] Bioconductor
# lattice                0.22-5    2023-10-24 [1] CRAN (R 4.3.1)
# limma                  3.58.1    2023-10-31 [1] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [1] CRAN (R 4.3.2)
# Matrix                 1.6-4     2023-11-30 [1] CRAN (R 4.3.2)
# MatrixGenerics       * 1.14.0    2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.2)
# R.methodsS3            1.8.2     2022-06-13 [1] CRAN (R 4.3.2)
# R.oo                   1.26.0    2024-01-24 [1] CRAN (R 4.3.2)
# R.utils                2.12.3    2023-11-18 [1] CRAN (R 4.3.2)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.2)
# RCurl                  1.98-1.13 2023-11-02 [1] CRAN (R 4.3.2)
# restfulr               0.0.15    2022-06-16 [1] CRAN (R 4.3.2)
# rhdf5                  2.46.1    2023-11-29 [1] Bioconductor 3.18 (R 4.3.2)
# rhdf5filters           1.14.1    2023-11-06 [1] Bioconductor
# Rhdf5lib               1.24.1    2023-12-11 [1] Bioconductor 3.18 (R 4.3.2)
# rjson                  0.2.21    2022-01-09 [1] CRAN (R 4.3.2)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.2)
# Rsamtools              2.18.0    2023-10-24 [1] Bioconductor
# rtracklayer          * 1.62.0    2023-10-24 [1] Bioconductor
# S4Arrays               1.2.0     2023-10-24 [1] Bioconductor
# S4Vectors            * 0.40.2    2023-11-23 [1] Bioconductor 3.18 (R 4.3.2)
# scuttle                1.12.0    2023-10-24 [1] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.2)
# SingleCellExperiment * 1.24.0    2023-10-24 [1] Bioconductor
# SparseArray            1.2.2     2023-11-07 [1] Bioconductor
# sparseMatrixStats      1.14.0    2023-10-24 [1] Bioconductor
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.2)
# SummarizedExperiment * 1.32.0    2023-10-24 [1] Bioconductor
# XML                    3.99-0.16 2023-11-29 [1] CRAN (R 4.3.2)
# XVector                0.42.0    2023-10-24 [1] Bioconductor
# yaml                   2.3.8     2023-12-11 [1] CRAN (R 4.3.2)
# zlibbioc               1.48.0    2023-10-24 [1] Bioconductor
# 
# [1] /jhpce/shared/libd/core/r_nac/1.0/nac_env/lib/R/library
# 
# ───────────────────────────────────────────────────────────────────────────────────────────────────────────
