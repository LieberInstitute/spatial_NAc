#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(dplyr)
library(scran)
library(here)

sce <- readRDS(here("processed-data","12_snRNA","sce_emptyDrops_removed_withQC.Rds"))

sce
# class: SingleCellExperiment 
# dim: 36601 120449 
# metadata(1): Samples
# assays(1): counts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(6): source type ... gene_name gene_type
# colnames(120449): 1_AAACCCAAGACCAACG-1 1_AAACCCACAGTCAGCC-1 ...
# 20_TTTGTTGCAAGATGTA-1 20_TTTGTTGGTACGAAAT-1
# colData names(30): Sample Barcode ... IEG_cells discard
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):

library(patchwork)

#By Sample 

Number_Genes_Sample <- plotColData(sce,x = "Sample",y = "detected",color_by = "discard") + 
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Number of Genes")

ggsave(plot = Number_Genes_Sample, 
       filename = here("plots","12_snRNA","Supplementary","QC_Supp_Figure","Numer_Genes_Sample.pdf"),
       height = 8,width = 8)

Library_Size_Sample <- plotColData(sce,x = "Sample",y = "sum",color_by = "discard") + 
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Library Size")

ggsave(plot = Library_Size_Sample, 
       filename = here("plots","12_snRNA","Supplementary","QC_Supp_Figure","Library_Size_Sample.pdf"),
       height = 8,width = 8)

ggsave(filename = here("plots","12_snRNA","Supplementary","QC_Supp_Figure","Numer_Genes_Sample.pdf"),height = 8,width = 8)

Percent_mito_Sample <- plotColData(sce,x = "Sample",y = "subsets_Mito_percent",color_by = "discard") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 5) +
  labs(y = "% Mitochondria")


ggsave(plot = Percent_mito_Sample, 
       filename = here("plots","12_snRNA","Supplementary","QC_Supp_Figure","Percent_mito_Sample.pdf"),
       height = 8,width = 8)


doublet_score_Sample <- plotColData(sce,x = "Sample",y = "doubletScore",color_by = "discard") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 5) 

ggsave(plot = doublet_score_Sample, 
       filename = here("plots","12_snRNA","Supplementary","QC_Supp_Figure","doublet_score_Sample.pdf"),
       height = 8,width = 8)


Sample_Figure_Only <- (Number_Genes | Library_Size) / ( Percent_mito| doublet_score ) + plot_annotation(tag_levels = "A")
ggsave(Sample_Figure_Only,
       filename = here("plots","12_snRNA","Supplementary","QC_Supp_Figure","Supp_Figure_SamplesOnly.pdf"),
       width = 10, 
       height = 10)



sessioninfo::session_info()
# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2025-07-03
# pandoc   3.1.3 @ /jhpce/shared/libd/core/r_nac/1.0/nac_env/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version    date (UTC) lib source
# abind                  1.4-5      2016-07-21 [1] CRAN (R 4.3.2)
# beachmat               2.18.0     2023-10-24 [1] Bioconductor
# beeswarm               0.4.0      2021-06-01 [1] CRAN (R 4.3.2)
# Biobase              * 2.62.0     2023-10-24 [1] Bioconductor
# BiocGenerics         * 0.48.1     2023-11-01 [1] Bioconductor
# BiocNeighbors          1.20.2     2024-01-07 [1] Bioconductor 3.18 (R 4.3.2)
# BiocParallel           1.36.0     2023-10-24 [1] Bioconductor
# BiocSingular           1.18.0     2023-10-24 [1] Bioconductor
# bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.3.2)
# bluster                1.11.4     2024-02-02 [1] Github (LTLA/bluster@17dd9c8)
# cli                    3.6.2      2023-12-11 [1] CRAN (R 4.3.2)
# cluster                2.1.6      2023-12-01 [1] CRAN (R 4.3.2)
# codetools              0.2-19     2023-02-01 [1] CRAN (R 4.3.0)
# colorspace             2.1-0      2023-01-23 [1] CRAN (R 4.3.0)
# cowplot                1.1.1      2020-12-30 [1] CRAN (R 4.3.2)
# crayon                 1.5.2      2022-09-29 [1] CRAN (R 4.3.0)
# DelayedArray           0.28.0     2023-10-24 [1] Bioconductor
# DelayedMatrixStats     1.24.0     2023-10-24 [1] Bioconductor
# dplyr                * 1.1.4      2023-11-17 [1] CRAN (R 4.3.2)
# dqrng                  0.3.2      2023-11-29 [1] CRAN (R 4.3.2)
# edgeR                  4.0.3      2023-12-10 [1] Bioconductor 3.18 (R 4.3.2)
# fansi                  1.0.6      2023-12-08 [1] CRAN (R 4.3.2)
# farver                 2.1.1      2022-07-06 [1] CRAN (R 4.3.0)
# generics               0.1.3      2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.38.1     2023-11-08 [1] Bioconductor
# GenomeInfoDbData       1.2.11     2023-12-12 [1] Bioconductor
# GenomicRanges        * 1.54.1     2023-10-29 [1] Bioconductor
# ggbeeswarm             0.7.2      2023-04-29 [1] CRAN (R 4.3.2)
# ggplot2              * 3.5.1      2024-04-23 [1] CRAN (R 4.3.2)
# ggrepel                0.9.4      2023-10-13 [1] CRAN (R 4.3.2)
# glue                   1.7.0      2024-01-09 [1] CRAN (R 4.3.2)
# gridExtra              2.3        2017-09-09 [1] CRAN (R 4.3.2)
# gtable                 0.3.4      2023-08-21 [1] CRAN (R 4.3.1)
# here                 * 1.0.1      2020-12-13 [1] CRAN (R 4.3.2)
# igraph                 2.0.3      2024-03-13 [1] CRAN (R 4.3.2)
# IRanges              * 2.36.0     2023-10-24 [1] Bioconductor
# irlba                  2.3.5.1    2022-10-03 [1] CRAN (R 4.3.2)
# labeling               0.4.3      2023-08-29 [1] CRAN (R 4.3.1)
# lattice                0.22-5     2023-10-24 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4      2023-11-07 [1] CRAN (R 4.3.2)
# limma                  3.58.1     2023-10-31 [1] Bioconductor
# locfit                 1.5-9.8    2023-06-11 [1] CRAN (R 4.3.2)
# magrittr               2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-4      2023-11-30 [1] CRAN (R 4.3.2)
# MatrixGenerics       * 1.14.0     2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0      2023-12-11 [1] CRAN (R 4.3.2)
# metapod                1.10.0     2023-10-24 [1] Bioconductor
# munsell                0.5.0      2018-06-12 [1] CRAN (R 4.3.0)
# patchwork            * 1.3.0.9000 2024-11-06 [1] Github (thomasp85/patchwork@2695a9f)
# pillar                 1.9.0      2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
# R6                     2.5.1      2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.7      2023-12-11 [1] CRAN (R 4.3.2)
# Rcpp                   1.0.12     2024-01-09 [1] CRAN (R 4.3.2)
# RCurl                  1.98-1.13  2023-11-02 [1] CRAN (R 4.3.2)
# rlang                  1.1.3      2024-01-10 [1] CRAN (R 4.3.2)
# rprojroot              2.0.4      2023-11-05 [1] CRAN (R 4.3.2)
# rsvd                   1.0.5      2021-04-16 [1] CRAN (R 4.3.2)
# S4Arrays               1.2.0      2023-10-24 [1] Bioconductor
# S4Vectors            * 0.40.2     2023-11-23 [1] Bioconductor 3.18 (R 4.3.2)
# ScaledMatrix           1.10.0     2023-10-24 [1] Bioconductor
# scales                 1.3.0      2023-11-28 [1] CRAN (R 4.3.2)
# scater               * 1.30.1     2023-11-16 [1] Bioconductor
# scran                * 1.30.0     2023-10-24 [1] Bioconductor
# scuttle              * 1.12.0     2023-10-24 [1] Bioconductor
# sessioninfo          * 1.2.2      2021-12-06 [1] CRAN (R 4.3.2)
# SingleCellExperiment * 1.24.0     2023-10-24 [1] Bioconductor
# SparseArray            1.2.2      2023-11-07 [1] Bioconductor
# sparseMatrixStats      1.14.0     2023-10-24 [1] Bioconductor
# statmod                1.5.0      2023-01-06 [1] CRAN (R 4.3.2)
# SummarizedExperiment * 1.32.0     2023-10-24 [1] Bioconductor
# systemfonts            1.2.1      2025-01-20 [1] CRAN (R 4.3.2)
# textshaping            0.3.7      2023-10-09 [1] CRAN (R 4.3.1)
# tibble                 3.2.1      2023-03-20 [1] CRAN (R 4.3.0)
# tidyselect             1.2.0      2022-10-10 [1] CRAN (R 4.3.0)
# utf8                   1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5      2023-12-01 [1] CRAN (R 4.3.2)
# vipor                  0.4.5      2017-03-22 [1] CRAN (R 4.3.2)
# viridis                0.6.4      2023-07-22 [1] CRAN (R 4.3.2)
# viridisLite            0.4.2      2023-05-02 [1] CRAN (R 4.3.0)
# withr                  2.5.2      2023-10-30 [1] CRAN (R 4.3.1)
# XVector                0.42.0     2023-10-24 [1] Bioconductor
# zlibbioc               1.48.0     2023-10-24 [1] Bioconductor
# 
# [1] /jhpce/shared/libd/core/r_nac/1.0/nac_env/lib/R/library
# 
# ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
