library("SpatialExperiment")
library("HDF5Array")
library("lobstr")
library("here")
library("sessioninfo")

spe_path_in <- here(
    "processed-data", "05_harmony_BayesSpace", "spe_filtered.rds"
)
spe_dir_out <- here("code", "06_deploy_app", "spe_shiny")

spe <- readRDS(spe_path_in)
spe
# class: SpatialExperiment
# dim: 30854 135597
# metadata(0):
# assays(2): counts logcounts
# rownames(30854): ENSG00000243485 ENSG00000238009 ... ENSG00000278817
#   ENSG00000277196
# rowData names(7): source type ... gene_type gene_search
# colnames(135597): AAACAACGAATAGTTC-1_V11D01-384_B1
#   AAACAAGTATCTCCCA-1_V11D01-384_B1 ... TTGTTTGTATTACACG-1_V13M06-378_D1
#   TTGTTTGTGTAAATTC-1_V13M06-378_D1
# colData names(41): sample_id in_tissue ... scran_quick_cluster
#   sizeFactor
# reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor
obj_size(spe)
# 6.99 GB

## Check what images we have stored
imgData(spe)
# DataFrame with 11 rows and 4 columns
#      sample_id    image_id   data scaleFactor
#    <character> <character> <list>   <numeric>
# 1       Br2720      lowres   ####   0.0215889
# 2       Br8492      lowres   ####   0.0281657
# 3       Br6522      lowres   ####   0.0235248
# 4       Br6423      lowres   ####   0.0271168
# 5       Br8325      lowres   ####   0.0200040
# 6       Br6471      lowres   ####   0.0188691
# 7       Br2743      lowres   ####   0.0266217
# 8       Br6432      lowres   ####   0.0230150
# 9       Br6432       hires   ####   0.0767165
# 10      Br6432    detected   ####   0.0767165
# 11      Br6432     aligned   ####   0.0767165

## Drop images we won't use
imgData(spe) <- imgData(spe)[imgData(spe)$image_id == "lowres", ]
obj_size(spe)
# 6.88 GB

## Drop counts assay
assays(spe)$counts <- NULL

## Check size in memory
obj_size(spe)
# 3.52 GB
## That's still too large for shinyapps

## Let's find out how many spots express each gene
gene_expr <- rowSums(logcounts(spe) > 0)

## How many genes are expressed in 10 or more spots?
table(gene_expr > 10)
# FALSE  TRUE
#  5561 25293

## 20 spots?
table(gene_expr > 20)
# FALSE  TRUE
#  7053 23801

## 30 spots?
table(gene_expr > 30)
# FALSE  TRUE
#  8018 22836

## 100 spots?
table(gene_expr > 100)
# FALSE  TRUE
# 10794 20060

## 1000 spots?
table(gene_expr > 1000)
# FALSE  TRUE 
# 16018 14836

## 10000 spots?
table(gene_expr > 10000)
# FALSE  TRUE 
# 23085  7769 

obj_size(logcounts(spe))
# 3.38 GB
obj_size(logcounts(spe[gene_expr > 1000, ]))
# 3.35 GB
obj_size(logcounts(spe[gene_expr > 10000, ]))
# 2.97 GB
## We have to go to genes expressed in at least 10,000 spots to drop under 3 GB
## for the logcounts

## Let's try first without doing so
# spe <- spe[gene_expr > 10000, ]
obj_size(spe)
# 3.52 GB

spe
# class: SpatialExperiment 
# dim: 30854 135597 
# metadata(0):
# assays(1): logcounts
# rownames(30854): ENSG00000243485 ENSG00000238009 ... ENSG00000278817 ENSG00000277196
# rowData names(7): source type ... gene_type gene_search
# colnames(135597): AAACAACGAATAGTTC-1_V11D01-384_B1 AAACAAGTATCTCCCA-1_V11D01-384_B1 ...
#   TTGTTTGTATTACACG-1_V13M06-378_D1 TTGTTTGTGTAAATTC-1_V13M06-378_D1
# colData names(41): sample_id in_tissue ... scran_quick_cluster sizeFactor
# reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor
saveHDF5SummarizedExperiment(spe, spe_dir_out, replace = TRUE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.1 Patched (2023-09-05 r85073)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2023-09-05
#  pandoc   3.1.1 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.0)
#  beachmat               2.16.0    2023-04-25 [2] Bioconductor
#  Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
#  BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
#  BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
#  bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.0)
#  cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.0)
#  codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
#  colorout               1.2-2     2023-05-06 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.0)
#  crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.0)
#  DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
#  DelayedMatrixStats     1.22.6    2023-08-28 [2] Bioconductor
#  digest                 0.6.33    2023-07-07 [2] CRAN (R 4.3.1)
#  dplyr                  1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
#  dqrng                  0.3.1     2023-08-30 [2] CRAN (R 4.3.1)
#  DropletUtils           1.20.0    2023-04-25 [2] Bioconductor
#  edgeR                  3.42.4    2023-05-31 [2] Bioconductor
#  fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.0)
#  fastmap                1.1.1     2023-02-24 [2] CRAN (R 4.3.0)
#  generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.0)
#  GenomeInfoDb         * 1.36.2    2023-08-25 [2] Bioconductor
#  GenomeInfoDbData       1.2.10    2023-04-11 [2] Bioconductor
#  GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
#  ggplot2                3.4.3     2023-08-14 [2] CRAN (R 4.3.1)
#  glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.0)
#  gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
#  HDF5Array              1.28.1    2023-05-01 [2] Bioconductor
#  here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.0)
#  htmltools              0.5.6     2023-08-10 [2] CRAN (R 4.3.1)
#  htmlwidgets            1.6.2     2023-03-17 [2] CRAN (R 4.3.0)
#  httpuv                 1.6.11    2023-05-11 [2] CRAN (R 4.3.0)
#  IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
#  jsonlite               1.8.7     2023-06-29 [2] CRAN (R 4.3.1)
#  later                  1.3.1     2023-05-02 [2] CRAN (R 4.3.0)
#  lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
#  lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.0)
#  limma                  3.56.2    2023-06-04 [2] Bioconductor
#  lobstr               * 1.1.2     2022-06-22 [2] CRAN (R 4.3.0)
#  locfit                 1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
#  magick                 2.7.5     2023-08-07 [2] CRAN (R 4.3.1)
#  magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.0)
#  Matrix                 1.6-1     2023-08-14 [3] CRAN (R 4.3.1)
#  MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
#  matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.0)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.0)
#  pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.0)
#  pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.0)
#  png                    0.1-8     2022-11-29 [2] CRAN (R 4.3.0)
#  prettyunits            1.1.1     2020-01-24 [2] CRAN (R 4.3.0)
#  promises               1.2.1     2023-08-10 [2] CRAN (R 4.3.1)
#  R.methodsS3            1.8.2     2022-06-13 [2] CRAN (R 4.3.0)
#  R.oo                   1.25.0    2022-06-12 [2] CRAN (R 4.3.0)
#  R.utils                2.12.2    2022-11-11 [2] CRAN (R 4.3.0)
#  R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.0)
#  Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
#  RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.0)
#  rhdf5                  2.44.0    2023-04-25 [2] Bioconductor
#  rhdf5filters           1.12.1    2023-04-30 [2] Bioconductor
#  Rhdf5lib               1.22.0    2023-04-25 [2] Bioconductor
#  rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.3.0)
#  rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.0)
#  rmote                  0.3.4     2023-05-06 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.0)
#  S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
#  S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
#  scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.0)
#  scuttle                1.10.2    2023-08-03 [2] Bioconductor
#  servr                  0.27      2023-05-02 [1] CRAN (R 4.3.0)
#  sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.0)
#  SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
#  sparseMatrixStats      1.12.2    2023-07-02 [2] Bioconductor
#  SpatialExperiment    * 1.10.0    2023-04-25 [2] Bioconductor
#  SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
#  tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.0)
#  tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.0)
#  utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.0)
#  vctrs                  0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
#  xfun                   0.40      2023-08-09 [2] CRAN (R 4.3.1)
#  XVector                0.40.0    2023-04-25 [2] Bioconductor
#  zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
# 
#  [1] /users/lcollado/R/4.3
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/R/4.3/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/R/4.3/lib64/R/library
# 
# ───────────────────────────────────────────────────────────────────────────────────────────────────────────
