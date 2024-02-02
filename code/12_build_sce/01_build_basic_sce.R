#Concatenate samples and build basic sce
#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc

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
#Read 10x data and create sce - 2024-02-01 22:27:13.350234
sce <- read10xCounts(samples = sample_data$Raw_data_path,
                     sample.names = sample_data$Sample,
                     type = "sparse",
                     col.names = TRUE)
message("RDone - ", Sys.time())
#RDone - 2024-02-01 22:35:33.506057

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
     file = here("processed-data","sce_raw.rda"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2024-02-01 23:13:54 EST"
# user   system  elapsed 
# 1420.613   36.826 4680.702 
# ─ Session info ───────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-02-01
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
# beachmat               2.16.0    2023-04-25 [2] Bioconductor
# Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
# BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
# BiocIO                 1.10.0    2023-04-25 [2] Bioconductor
# BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
# Biostrings             2.68.1    2023-05-16 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
# colorout             * 1.3-0.1   2023-12-01 [1] Github (jalvesaq/colorout@deda341)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
# DelayedMatrixStats     1.22.6    2023-08-28 [2] Bioconductor
# dqrng                  0.3.1     2023-08-30 [2] CRAN (R 4.3.1)
# DropletUtils         * 1.20.0    2023-04-25 [2] Bioconductor
# edgeR                  3.42.4    2023-05-31 [2] Bioconductor
# GenomeInfoDb         * 1.36.3    2023-09-07 [2] Bioconductor
# GenomeInfoDbData       1.2.10    2023-07-20 [2] Bioconductor
# GenomicAlignments      1.36.0    2023-04-25 [2] Bioconductor
# GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
# HDF5Array              1.28.1    2023-05-01 [2] Bioconductor
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
# lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# limma                  3.56.2    2023-06-04 [2] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
# Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# R.methodsS3            1.8.2     2022-06-13 [2] CRAN (R 4.3.1)
# R.oo                   1.25.0    2022-06-12 [2] CRAN (R 4.3.1)
# R.utils                2.12.2    2022-11-11 [2] CRAN (R 4.3.1)
# Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.3.1)
# rhdf5                  2.44.0    2023-04-25 [2] Bioconductor
# rhdf5filters           1.12.1    2023-04-30 [2] Bioconductor
# Rhdf5lib               1.22.1    2023-09-10 [2] Bioconductor
# rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# Rsamtools              2.16.0    2023-04-25 [2] Bioconductor
# rtracklayer          * 1.60.1    2023-08-15 [2] Bioconductor
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# scuttle                1.10.2    2023-08-03 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
# sparseMatrixStats      1.12.2    2023-07-02 [2] Bioconductor
# SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
# XML                    3.99-0.14 2023-03-19 [2] CRAN (R 4.3.1)
# XVector                0.40.0    2023-04-25 [2] Bioconductor
# yaml                   2.3.7     2023-01-23 [2] CRAN (R 4.3.1)
# zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
# 
# [1] /users/rphillip/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────────
