#Completed in interactive job: srun --cpus-per-task=2 --mem=20G --pty --x11 bash
#cd cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc
#module load r_nac
#code modified from: https://github.com/LieberInstitute/spatialdACC/tree/main/code/snRNA-seq/06_NMF

library(SingleCellExperiment)
library(pheatmap)
library(reshape2)
library(here)

#Load NMF results and the sce object
nmf_res <- readRDS(here("processed-data","12_snRNA","NMF","NMF_Results.Rds"))
sce <-  readRDS(here("processed-data","12_snRNA","sce_CellType_noresiduals.Rds"))

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

###########Correlate Brain ID with NMF patterns#########
##onehot encoding of the brain variable##
#Make a dataframe of brain IDs 
data <- as.data.frame(sce$Brain_ID)
colnames(data) <- "Brain_ID"

head(data)
# Brain_ID
# 1   Br8325
# 2   Br8325
# 3   Br8325
# 4   Br8325
# 5   Br8325
# 6   Br8325

dim(data)
#[1] 103785      1

onehot_brain <-  dcast(data = data, rownames(data) ~ Brain_ID, length)

dim(onehot_brain)
#[1] 103785     11

head(onehot_brain)
# rownames(data) Br2720 Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325 Br8492
# 1              1      0      0      0      0      0      0      0      1      0
# 2             10      0      0      0      0      0      0      0      1      0
# 3            100      0      0      0      0      0      0      0      1      0
# 4           1000      0      0      0      0      0      0      0      1      0
# 5          10000      0      0      0      0      0      0      0      0      1
# 6         100000      0      0      0      0      0      0      0      0      0
# Br8667
# 1      0
# 2      0
# 3      0
# 4      0
# 5      0
# 6      1

#Make the rownames the first column
rownames(onehot_brain) <- onehot_brain[,1]

#Convert first column to numeric 
onehot_brain[,1]<-as.numeric(onehot_brain[,1])

#Reorder based on first column
onehot_brain <- onehot_brain[order(onehot_brain[,1],decreasing=FALSE),]

head(onehot_brain)
# rownames(data) Br2720 Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325 Br8492
# 1              1      0      0      0      0      0      0      0      1      0
# 2              2      0      0      0      0      0      0      0      1      0
# 3              3      0      0      0      0      0      0      0      1      0
# 4              4      0      0      0      0      0      0      0      1      0
# 5              5      0      0      0      0      0      0      0      1      0
# 6              6      0      0      0      0      0      0      0      1      0
# Br8667
# 1      0
# 2      0
# 3      0
# 4      0
# 5      0
# 6      0

#Remove the first column
onehot_brain[,1] <- NULL

#Correlate brain ID with NMF patterns. 
pdf(here("plots","12_snRNA","NMF", "NMF_BrainID_correlation_heatmap.pdf"))
pheatmap(cor(t(nmf_res@h),onehot_brain), fontsize_row = 5)
dev.off()

###########Correlate QC measures with NMF patterns#########
##onehot encoding of the brain variable##
#Make a dataframe of brain IDs 
QC_vars <- colData(sce)[,c("detected","sum","subsets_Mito_percent")]

head(QC_vars)
# DataFrame with 6 rows and 3 columns
# detected       sum subsets_Mito_percent
# <integer> <numeric>            <numeric>
# 1_AAACCCAAGACCAACG-1      1744      3673            0.0544514
# 1_AAACCCACAGTCAGCC-1      2440      5977            0.0000000
# 1_AAACCCATCACCCTCA-1      6355     21932            0.0775123
# 1_AAACCCATCTACGGGC-1      1845      3174            0.1260239
# 1_AAACGAAAGGTGAGCT-1      1918      3560            0.0280899
1_AAACGAACATAGACTC-1      6795     27502            0.0145444

#Make the detected numeric
QC_vars$detected <- as.numeric(QC_vars$detected)

#Convert to matrix
QC_vars <- as.matrix(QC_vars)

pdf(here("plots","12_snRNA","NMF","NMF_QC_correlation_heatmap.pdf"))
pheatmap(cor(t(nmf_res@h),QC_vars), fontsize_row = 5)
dev.off()


###########Correlate CellType with NMF patterns#########
#Create datafarme of celltype 
ct_data <- as.data.frame(sce$CellType.Final)
colnames(ct_data) <- "CellType.Final" 

#One hot encode the cell type
onehot_CellType <-  dcast(data = ct_data, rownames(ct_data) ~ CellType.Final, length)

head(onehot_CellType)
# rownames(ct_data) Astrocyte_A Astrocyte_B DRD1_MSN_A DRD1_MSN_B DRD1_MSN_C
# 1                 1           0           0          0          0          0
# 2                10           0           0          0          0          1
# 3               100           0           0          1          0          0
# 4              1000           0           0          0          0          0
# 5             10000           0           0          0          0          0
# 6            100000           0           0          0          0          0
# DRD1_MSN_D DRD2_MSN_A DRD2_MSN_B Ependymal Excitatory Inh_A Inh_B Inh_C Inh_D
# 1          0          0          0         0          0     0     0     0     0
# 2          0          0          0         0          0     0     0     0     0
# 3          0          0          0         0          0     0     0     0     0
# 4          0          0          1         0          0     0     0     0     0
# 5          0          0          1         0          0     0     0     0     0
# 6          0          0          1         0          0     0     0     0     0
# Inh_E Inh_F Inh_G Inh_H Inh_I Macrophage Microglia
# 1     0     0     0     0     0          0         0
# 2     0     0     0     0     0          0         0
# 3     0     0     0     0     0          0         0
# 4     0     0     0     0     0          0         0
# 5     0     0     0     0     0          0         0
# 6     0     0     0     0     0          0         0
# Mural_Endothelial_Fibroblast Oligo OPC T-Cell
# 1                            0     1   0      0
# 2                            0     0   0      0
# 3                            0     0   0      0
# 4                            0     0   0      0
# 5                            0     0   0      0
# 6                            0     0   0      0

rownames(onehot_CellType) <- onehot_CellType[,1]
onehot_CellType[,1] <- as.numeric(onehot_CellType[,1])
onehot_CellType <- onehot_CellType[order(onehot_CellType[,1],decreasing=FALSE),]

head(onehot_CellType)
# rownames(ct_data) Astrocyte_A Astrocyte_B DRD1_MSN_A DRD1_MSN_B DRD1_MSN_C
# 1                 1           0           0          0          0          0
# 2                 2           0           0          0          0          0
# 3                 3           0           0          0          0          1
# 4                 4           0           0          0          0          0
# 5                 5           0           0          0          0          0
# 6                 6           0           0          0          0          1
# DRD1_MSN_D DRD2_MSN_A DRD2_MSN_B Ependymal Excitatory Inh_A Inh_B Inh_C Inh_D
# 1          0          0          0         0          0     0     0     0     0
# 2          0          0          0         0          0     0     0     0     0
# 3          0          0          0         0          0     0     0     0     0
# 4          0          0          0         0          0     0     0     0     0
# 5          0          0          0         0          0     0     0     0     0
# 6          0          0          0         0          0     0     0     0     0
# Inh_E Inh_F Inh_G Inh_H Inh_I Macrophage Microglia
# 1     0     0     0     0     0          0         0
# 2     0     0     0     0     0          0         0
# 3     0     0     0     0     0          0         0
# 4     0     0     0     0     0          0         1
# 5     0     0     0     0     0          1         0
# 6     0     0     0     0     0          0         0
# Mural_Endothelial_Fibroblast Oligo OPC T-Cell
# 1                            0     1   0      0
# 2                            0     1   0      0
# 3                            0     0   0      0
# 4                            0     0   0      0
# 5                            0     0   0      0
# 6                            0     0   0      0

onehot_CellType[,1]<-NULL


###correlate with nmf patterns
pdf(here("plots","12_snRNA","NMF","NMF_CellType_correlation_heatmap.pdf"))
pheatmap(cor(t(nmf_res@h),onehot_CellType), fontsize_row = 5)
dev.off()

########### Aggregate NMF patterns #########
# create dataframe
aggr_data <- data.frame(colData(sce), t(nmf_res@h))

# aggregate NMF patterns across cell types
# grep "NMF" to get all NMF patterns
aggr_data2 <- aggregate(x = aggr_data[,grep("nmf", colnames(aggr_data))],
                        by = list(aggr_data$CellType.Final),
                        FUN = mean)

aggr_data2[1:5,1:5]
# Group.1         nmf1         nmf2         nmf3         nmf4
# 1 Astrocyte_A 1.202200e-05 5.851288e-08 8.490542e-08 2.945811e-06
# 2 Astrocyte_B 9.471367e-06 4.281921e-08 5.002034e-08 2.788676e-06
# 3  DRD1_MSN_A 1.200896e-05 5.516762e-07 9.030000e-07 1.318676e-05
# 4  DRD1_MSN_B 1.045089e-05 2.072535e-06 1.241188e-05 9.116840e-06
# 5  DRD1_MSN_C 7.855158e-06 2.406406e-05 3.603876e-05 1.644265e-05
#aggr_data2 contains mean nmf values by CellType.Final

# move Group.1 to row names, then drop
rownames(aggr_data2) <- aggr_data2$Group.1
aggr_data2 <- aggr_data2[,-1]

pdf(here("plots","12_snRNA","NMF","NMF_CellType_correlation_aggregated_heatmap.pdf"))
pheatmap(aggr_data2,
               color=colorRampPalette(c("blue","white","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               scale="column",
               fontsize_col = 5)
dev.off()



#Pull top 50 genes for each NMF pattern
#Make w (feature factor matrix) a separate dataframe
w_df <- as.data.frame(nmf_res@w)

#Make an ensemble column and then merge with gene name info
w_df$ensembl <- rownames(w_df)
w_df_merge <- merge(x = w_df,
                    y = rowData(sce)[,c("gene_id","gene_name")],
                    by.x = "ensembl",
                    by.y = "gene_id")

#Change rownames to the gene_name
rownames(w_df_merge) <- w_df_merge$gene_name

#Keep only nmf patterns
w_df_merge <- as.matrix(w_df_merge[,c(2:67)])
 
# function for getting top n genes for each pattern
top_genes <- function(W, n=50){
  top_genes <- apply(W, 2, function(x) names(sort(x, decreasing=TRUE)[1:n]))
  return(top_genes)
}

# get top 10 genes
top50 <- top_genes(w_df_merge, 50)
write.csv(top50, 
          file = here("processed-data","12_snRNA","NMF", "NMF_top50_genes.csv"),
          quote = FALSE)

#lots of the top 10 features in each nmf pattern. 
library(scater)
library(ggplot2)

load(here("processed-data","12_snRNA","CellType_Final_Cluster_Cols.rda"),verbose = TRUE)
# Loading objects:
#   cluster_cols

#just to see, check top 10 genes for nmf27 which corresponds to one of the D1-island clusters
plotExpression(object = sce,features = top50[1:10,"nmf27"],
               x = "CellType.Final",swap_rownames = "gene_name",
               color_by = "CellType.Final",ncol = 2) +
  scale_color_manual(values = cluster_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

##Add nmf patterns to sce object
loads <- t(nmf_res@h)

#Check that they are in the same order
all(rownames(loads) == rownames(colData(sce)))
#[1] TRUE

#Add the patterns 
colData(sce) <- cbind(colData(sce),loads)

for(i in 1:66){
  print(i)
  x <- ggplot(data = as.data.frame(colData(sce)),
              aes_string(x = "CellType.Final",
                         y=paste0("nmf",i),
                         fill = "CellType.Final")) +
    geom_boxplot() +
    scale_fill_manual(values = cluster_cols) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") 
  ggsave(plot = x,
         filename = here("plots","12_snRNA","NMF","nmf_celltype_vlns",
                         paste0("nmf",i,".png")),
         height = 6, 
         width = 8)
}

#session info
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
# [1] "Reproducibility information:"
# [1] "2024-05-21 17:14:19 EDT"
# user    system   elapsed 
# 3917.374    15.946 14740.936 
# ─ Session info ─────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-05-21
# pandoc   3.1.3 @ /jhpce/shared/libd/core/r_nac/1.0/nac_env/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────
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
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.2)
# codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.0)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# cowplot                1.1.1     2020-12-30 [1] CRAN (R 4.3.2)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# DelayedArray           0.28.0    2023-10-24 [1] Bioconductor
# DelayedMatrixStats     1.24.0    2023-10-24 [1] Bioconductor
# dplyr                  1.1.4     2023-11-17 [1] CRAN (R 4.3.2)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.2)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.38.1    2023-11-08 [1] Bioconductor
# GenomeInfoDbData       1.2.11    2023-12-12 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-29 [1] Bioconductor
# ggbeeswarm             0.7.2     2023-04-29 [1] CRAN (R 4.3.2)
# ggplot2              * 3.4.4     2023-10-12 [1] CRAN (R 4.3.2)
# ggrepel                0.9.4     2023-10-13 [1] CRAN (R 4.3.2)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.2)
# gridExtra              2.3       2017-09-09 [1] CRAN (R 4.3.2)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.1)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.2)
# IRanges              * 2.36.0    2023-10-24 [1] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.3.2)
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.1)
# lattice                0.22-5    2023-10-24 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.2)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-4     2023-11-30 [1] CRAN (R 4.3.2)
# MatrixGenerics       * 1.14.0    2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.2)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# pheatmap             * 1.0.12    2019-01-04 [1] CRAN (R 4.3.2)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# plyr                   1.8.9     2023-10-02 [1] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.7     2023-12-11 [1] CRAN (R 4.3.2)
# RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.2)
# RcppML               * 0.3.7     2021-09-21 [1] CRAN (R 4.3.2)
# RCurl                  1.98-1.13 2023-11-02 [1] CRAN (R 4.3.2)
# reshape2             * 1.4.4     2020-04-09 [1] CRAN (R 4.3.0)
# rlang                  1.1.3     2024-01-10 [1] CRAN (R 4.3.2)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.2)
# rsvd                   1.0.5     2021-04-16 [1] CRAN (R 4.3.2)
# S4Arrays               1.2.0     2023-10-24 [1] Bioconductor
# S4Vectors            * 0.40.2    2023-11-23 [1] Bioconductor 3.18 (R 4.3.2)
# ScaledMatrix           1.10.0    2023-10-24 [1] Bioconductor
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.2)
# scater               * 1.30.1    2023-11-16 [1] Bioconductor
# scuttle              * 1.12.0    2023-10-24 [1] Bioconductor
# sessioninfo            1.2.2     2021-12-06 [1] CRAN (R 4.3.2)
# SingleCellExperiment * 1.24.0    2023-10-24 [1] Bioconductor
# SparseArray            1.2.2     2023-11-07 [1] Bioconductor
# sparseMatrixStats      1.14.0    2023-10-24 [1] Bioconductor
# stringi                1.8.3     2023-12-11 [1] CRAN (R 4.3.2)
# stringr                1.5.1     2023-11-14 [1] CRAN (R 4.3.2)
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
# ────────────────────────────────────────────────────────────────────────────────────────
# 
