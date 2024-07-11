#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

library(SingleCellExperiment)
library(ComplexHeatmap)
library(sessioninfo)
library(ggplot2)
library(scater)
library(dplyr)
library(scran)
library(here)

sce <- readRDS(here("processed-data","12_snRNA","sce_CellType_noresiduals.Rds"))

sce
# class: SingleCellExperiment 
# dim: 36601 103785 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(103785): 1_AAACCCAAGACCAACG-1 1_AAACCCACAGTCAGCC-1 ...
# 20_TTTGTTGCAAGATGTA-1 20_TTTGTTGGTACGAAAT-1
# colData names(33): Sample Barcode ... sizeFactor CellType.Final
# reducedDimNames(5): GLMPCA_approx tSNE HARMONY tSNE_HARMONY
# umap_HARMONY
# mainExpName: NULL
# altExpNames(0):


#DRD1 featureplot
d1_Feature <- plotReducedDim(object = sce,
                             dimred = "umap_HARMONY",
                             color_by = "DRD1",
                             swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","red")) +
  theme(axis.line = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) 

#with legend
ggsave(plot = d1_Feature,
       filename = here("plots","12_snRNA","D1_FeaturePlot_legend.pdf"),
       height = 8,
       width = 8)
#without legend
ggsave(plot = d1_Feature + theme(legend.position = "none"),
       filename = here("plots","12_snRNA","D1_FeaturePlot_nolegend.png"),
       height = 8,
       width = 8)

#DRD2 featureplot
d2_Feature <- plotReducedDim(object = sce,
                             dimred = "umap_HARMONY",
                             color_by = "DRD2",
                             swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","red")) +
  theme(axis.line = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) 

#with legend
ggsave(plot = d2_Feature,
       filename = here("plots","12_snRNA","D2_FeaturePlot_legend.pdf"),
       height = 8,
       width = 8)
#without legend
ggsave(plot = d2_Feature + theme(legend.position = "none"),
       filename = here("plots","12_snRNA","D2_FeaturePlot_nolegend.png"),
       height = 8,
       width = 8)


print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
# [1] "Reproducibility information:"
# [1] "2024-07-11 13:18:07 EDT"
# user   system  elapsed 
# 117.278   13.433 6863.154 
# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-07-11
# pandoc   3.1.3 @ /jhpce/shared/libd/core/r_nac/1.0/nac_env/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.2)
# beachmat               2.18.0    2023-10-24 [1] Bioconductor
# beeswarm               0.4.0     2021-06-01 [1] CRAN (R 4.3.2)
# Biobase              * 2.62.0    2023-10-24 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-01 [1] Bioconductor
# BiocNeighbors          1.20.2    2024-01-07 [1] Bioconductor 3.18 (R 4.3.2)
# BiocParallel           1.36.0    2023-10-24 [1] Bioconductor
# BiocSingular           1.18.0    2023-10-24 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.2)
# bluster                1.11.4    2024-02-02 [1] Github (LTLA/bluster@17dd9c8)
# circlize               0.4.16    2024-02-20 [1] CRAN (R 4.3.2)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.2)
# clue                   0.3-65    2023-09-23 [1] CRAN (R 4.3.2)
# cluster                2.1.6     2023-12-01 [1] CRAN (R 4.3.2)
# codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.0)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# ComplexHeatmap       * 2.18.0    2023-10-24 [1] Bioconductor
# cowplot                1.1.1     2020-12-30 [1] CRAN (R 4.3.2)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# DelayedArray           0.28.0    2023-10-24 [1] Bioconductor
# DelayedMatrixStats     1.24.0    2023-10-24 [1] Bioconductor
# digest                 0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
# doParallel             1.0.17    2022-02-07 [1] CRAN (R 4.3.2)
# dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.3.2)
# dqrng                  0.3.2     2023-11-29 [1] CRAN (R 4.3.2)
# edgeR                  4.0.3     2023-12-10 [1] Bioconductor 3.18 (R 4.3.2)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.2)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# foreach                1.5.2     2022-02-02 [1] CRAN (R 4.3.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.38.1    2023-11-08 [1] Bioconductor
# GenomeInfoDbData       1.2.11    2023-12-12 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-29 [1] Bioconductor
# GetoptLong             1.0.5     2020-12-15 [1] CRAN (R 4.3.2)
# ggbeeswarm             0.7.2     2023-04-29 [1] CRAN (R 4.3.2)
# ggplot2              * 3.4.4     2023-10-12 [1] CRAN (R 4.3.2)
# ggrepel                0.9.4     2023-10-13 [1] CRAN (R 4.3.2)
# GlobalOptions          0.1.2     2020-06-10 [1] CRAN (R 4.3.2)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.2)
# gridExtra              2.3       2017-09-09 [1] CRAN (R 4.3.2)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.1)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.2)
# igraph                 2.0.3     2024-03-13 [1] CRAN (R 4.3.2)
# IRanges              * 2.36.0    2023-10-24 [1] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.3.2)
# iterators              1.0.14    2022-02-05 [1] CRAN (R 4.3.0)
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.1)
# lattice                0.22-5    2023-10-24 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.2)
# limma                  3.58.1    2023-10-31 [1] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [1] CRAN (R 4.3.2)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-4     2023-11-30 [1] CRAN (R 4.3.2)
# MatrixGenerics       * 1.14.0    2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.2)
# metapod                1.10.0    2023-10-24 [1] Bioconductor
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# png                    0.1-8     2022-11-29 [1] CRAN (R 4.3.2)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.7     2023-12-11 [1] CRAN (R 4.3.2)
# RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.2)
# RCurl                  1.98-1.13 2023-11-02 [1] CRAN (R 4.3.2)
# rjson                  0.2.21    2022-01-09 [1] CRAN (R 4.3.2)
# rlang                  1.1.3     2024-01-10 [1] CRAN (R 4.3.2)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.2)
# rsvd                   1.0.5     2021-04-16 [1] CRAN (R 4.3.2)
# S4Arrays               1.2.0     2023-10-24 [1] Bioconductor
# S4Vectors            * 0.40.2    2023-11-23 [1] Bioconductor 3.18 (R 4.3.2)
# ScaledMatrix           1.10.0    2023-10-24 [1] Bioconductor
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.2)
# scater               * 1.30.1    2023-11-16 [1] Bioconductor
# scran                * 1.30.0    2023-10-24 [1] Bioconductor
# scuttle              * 1.12.0    2023-10-24 [1] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.2)
# shape                  1.4.6     2021-05-19 [1] CRAN (R 4.3.0)
# SingleCellExperiment * 1.24.0    2023-10-24 [1] Bioconductor
# SparseArray            1.2.2     2023-11-07 [1] Bioconductor
# sparseMatrixStats      1.14.0    2023-10-24 [1] Bioconductor
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.2)
# SummarizedExperiment * 1.32.0    2023-10-24 [1] Bioconductor
# systemfonts            1.0.5     2023-10-09 [1] CRAN (R 4.3.1)
# textshaping            0.3.7     2023-10-09 [1] CRAN (R 4.3.1)
# tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.2)
# vipor                  0.4.5     2017-03-22 [1] CRAN (R 4.3.2)
# viridis                0.6.4     2023-07-22 [1] CRAN (R 4.3.2)
# viridisLite            0.4.2     2023-05-02 [1] CRAN (R 4.3.0)
# withr                  2.5.2     2023-10-30 [1] CRAN (R 4.3.1)
# XVector                0.42.0    2023-10-24 [1] Bioconductor
# zlibbioc               1.48.0    2023-10-24 [1] Bioconductor
# 
# [1] /jhpce/shared/libd/core/r_nac/1.0/nac_env/lib/R/library
# 
# ────────────────────────────────────────────────────────────────────────────────────────────────────
