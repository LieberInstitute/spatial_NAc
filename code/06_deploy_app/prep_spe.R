library("SpatialExperiment")
library("HDF5Array")
library("lobstr")
library("here")
library("sessioninfo")
library("Matrix")

spe_dir_in <- here(
    "processed-data", "05_harmony_BayesSpace", "spe_filtered_hdf5"
)
spe_dir_out <- here("code", "06_deploy_app", "spe_shiny")

spe <- loadHDF5SummarizedExperiment(spe_dir_in)
spe
# class: SpatialExperiment
# dim: 31333 177804
# metadata(0):
# assays(3): counts logcounts binomial_deviance_residuals
# rownames(31333): ENSG00000243485 ENSG00000238009 ... ENSG00000278817
#   ENSG00000277196
# rowData names(8): source type ... gene_search binomial_deviance
# colnames(177804): AAACAACGAATAGTTC-1_V11D01-384_A1
#   AAACAAGTATCTCCCA-1_V11D01-384_A1 ... TTGTTTGTATTACACG-1_V13M13-362_C1
#   TTGTTTGTGTAAATTC-1_V13M13-362_C1
# colData names(41): sample_id in_tissue ... overlap_slide
#   exclude_overlapping
# reducedDimNames(5): 10x_pca 10x_tsne 10x_umap PCA GLMPCA_approx
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

#   All images should be low-res (if others existed, we'd drop them)
stopifnot(all(imgData(spe)$image_id == "lowres"))

#   Only keep logcounts
assays(spe) <- list(logcounts = assays(spe)$logcounts)

obj_size(spe)
# 336.46 MB
# (because logcounts are on disk, but Shiny appears to respect on-disk operations)

# (Comment out original gene-filtering steps, as the HDF5-backed SPE seems to
# behave well on Shiny, not loading all gene counts into memory)

# ## Check size in memory
# obj_size(spe)
# # 3.52 GB
# ## That's still too large for shinyapps

# ## Let's find out how many spots express each gene
# gene_expr <- rowSums(logcounts(spe) > 0)

# ## How many genes are expressed in 10 or more spots?
# table(gene_expr > 10)
# # FALSE  TRUE
# #  5561 25293

# ## 20 spots?
# table(gene_expr > 20)
# # FALSE  TRUE
# #  7053 23801

# ## 30 spots?
# table(gene_expr > 30)
# # FALSE  TRUE
# #  8018 22836

# ## 100 spots?
# table(gene_expr > 100)
# # FALSE  TRUE
# # 10794 20060

# ## 1000 spots?
# table(gene_expr > 1000)
# # FALSE  TRUE
# # 16018 14836

# ## 10000 spots?
# table(gene_expr > 10000)
# # FALSE  TRUE
# # 23085  7769

# obj_size(logcounts(spe))
# # 3.38 GB
# obj_size(logcounts(spe[gene_expr > 1000, ]))
# # 3.35 GB
# obj_size(logcounts(spe[gene_expr > 10000, ]))
# # 2.97 GB
# ## We have to go to genes expressed in at least 10,000 spots to drop under 3 GB
# ## for the logcounts

# ## Let's try first without doing so
# # spe <- spe[gene_expr > 10000, ]
# obj_size(spe)
# # 3.52 GB

# spe
# # class: SpatialExperiment
# # dim: 30854 135597
# # metadata(0):
# # assays(1): logcounts
# # rownames(30854): ENSG00000243485 ENSG00000238009 ... ENSG00000278817 ENSG00000277196
# # rowData names(7): source type ... gene_type gene_search
# # colnames(135597): AAACAACGAATAGTTC-1_V11D01-384_B1 AAACAAGTATCTCCCA-1_V11D01-384_B1 ...
# #   TTGTTTGTATTACACG-1_V13M06-378_D1 TTGTTTGTGTAAATTC-1_V13M06-378_D1
# # colData names(41): sample_id in_tissue ... scran_quick_cluster sizeFactor
# # reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# # mainExpName: NULL
# # altExpNames(0):
# # spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# # imgData names(4): sample_id image_id data scaleFactor
saveHDF5SummarizedExperiment(spe, spe_dir_out, replace = TRUE)

## Reproducibility information
print("Reproducibility information:")
Sys.time() # "2023-12-01 10:13:02 EST"
proc.time()
#     user   system  elapsed
# 1568.030   13.879 5401.584
options(width = 120)
sessioninfo::session_info()
# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.1 Patched (2023-07-19 r84711)
#  os       Rocky Linux 9.2 (Blue Onyx)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2023-12-01
#  pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc

# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  abind                * 1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
#  beachmat               2.16.0    2023-04-25 [2] Bioconductor
#  Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
#  BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
#  BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
#  bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
#  cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
#  codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
#  crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
#  DelayedArray         * 0.26.7    2023-07-28 [2] Bioconductor
#  DelayedMatrixStats     1.22.6    2023-08-28 [2] Bioconductor
#  dqrng                  0.3.1     2023-08-30 [2] CRAN (R 4.3.1)
#  DropletUtils           1.20.0    2023-04-25 [2] Bioconductor
#  edgeR                  3.42.4    2023-05-31 [2] Bioconductor
#  GenomeInfoDb         * 1.36.3    2023-09-07 [2] Bioconductor
#  GenomeInfoDbData       1.2.10    2023-07-20 [2] Bioconductor
#  GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
#  HDF5Array            * 1.28.1    2023-05-01 [2] Bioconductor
#  here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
#  IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
#  lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
#  limma                  3.56.2    2023-06-04 [2] Bioconductor
#  lobstr               * 1.1.2     2022-06-22 [2] CRAN (R 4.3.1)
#  locfit                 1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
#  magick                 2.7.5     2023-08-07 [2] CRAN (R 4.3.1)
#  magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
#  Matrix               * 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
#  MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
#  matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
#  prettyunits            1.2.0     2023-09-24 [1] CRAN (R 4.3.1)
#  R.methodsS3            1.8.2     2022-06-13 [2] CRAN (R 4.3.1)
#  R.oo                   1.25.0    2022-06-12 [2] CRAN (R 4.3.1)
#  R.utils                2.12.2    2022-11-11 [2] CRAN (R 4.3.1)
#  Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
#  RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
#  rhdf5                * 2.44.0    2023-04-25 [2] Bioconductor
#  rhdf5filters           1.12.1    2023-04-30 [2] Bioconductor
#  Rhdf5lib               1.22.1    2023-09-10 [2] Bioconductor
#  rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.3.1)
#  rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
#  rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
#  S4Arrays             * 1.2.0     2023-10-24 [1] Bioconductor
#  S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
#  scuttle                1.10.2    2023-08-03 [2] Bioconductor
#  sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
#  SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
#  sparseMatrixStats      1.12.2    2023-07-02 [2] Bioconductor
#  SpatialExperiment    * 1.10.0    2023-04-25 [2] Bioconductor
#  SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
#  XVector                0.40.0    2023-04-25 [2] Bioconductor
#  zlibbioc               1.46.0    2023-04-25 [2] Bioconductor

#  [1] /users/neagles/R/4.3
#  [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
#  [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library

# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
