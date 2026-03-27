# cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
# module load r_nac
library(SingleCellExperiment)
library(ggplot2)
library(scater)
library(scran)
library(lisi)
library(here)


#Load the object 
#load the sce object
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
# colData names(41): Sample Barcode ... sizeFactor CellType.Final
# reducedDimNames(4): GLMPCA_approx tSNE HARMONY tSNE_HARMONY
# mainExpName: NULL
# altExpNames(0):


#Subset for Inh_D 
sce_sub <- sce[,sce$CellType.Final == "Inh_D"]

#Pull counts for genes
# ENSG00000115665  - SLC5A7
# ENSG00000070748 - CHAT
# ENSG00000128683 - GAD1 
# ENSG00000136750 - GAD2 
# ENSG00000104888 - SLC17A7
counts_mat <- counts(sce_sub)[c("ENSG00000115665","ENSG00000070748",
                                "ENSG00000128683","ENSG00000136750",
                                "ENSG00000104888"),]


pct_mat <- matrix(nrow = 5, ncol = 5)
rownames(pct_mat) <- rownames(counts_mat)
colnames(pct_mat) <- rownames(counts_mat)


for(i in colnames(pct_mat)){
  print(i)
  for(l in rownames(pct_mat)){
    print(l)
    if(i ==l){
      next
    }else{
      pct_mat[l,i] <- sum(counts_mat[i,]!=0 & counts_mat[l,]!=0)/ncol(counts_mat) * 100
    }
  }
}


pct_mat
# ENSG00000115665 ENSG00000070748 ENSG00000128683 ENSG00000136750
# ENSG00000115665              NA       85.398230       89.823009      34.5132743
# ENSG00000070748       85.398230              NA       82.300885      30.9734513
# ENSG00000128683       89.823009       82.300885              NA      34.0707965
# ENSG00000136750       34.513274       30.973451       34.070796              NA
# ENSG00000104888        2.654867        2.654867        2.654867       0.4424779
# ENSG00000104888
# ENSG00000115665       2.6548673
# ENSG00000070748       2.6548673
# ENSG00000128683       2.6548673
# ENSG00000136750       0.4424779
# ENSG00000104888              NA

colnames(pct_mat) <- c("SLC5A7","CHAT","GAD1","GAD2","SLC17A7")
rownames(pct_mat) <- colnames(pct_mat)

gt <-pheatmap::pheatmap(
  pct_mat,
  na_col = "gray",
  cluster_rows = FALSE,
  cluster_cols = FALSE
)$gtable
ggsave(here("plots","12_snRNA","Inh_D_coexpression.pdf"), plot=gt)



###Reproduciblity
sessionInfo()
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-conda-linux-gnu (64-bit)
# Running under: Rocky Linux 9.4 (Blue Onyx)
# 
# Matrix products: default
# BLAS/LAPACK: /jhpce/shared/libd/core/r_nac/1.0/nac_env/lib/libopenblasp-r0.3.25.so;  LAPACK version 3.11.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: US/Eastern
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
#   [1] here_1.0.1                  lisi_1.0                   
# [3] scran_1.30.0                scater_1.30.1              
# [5] scuttle_1.12.0              ggplot2_3.5.1              
# [7] SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
# [9] Biobase_2.62.0              GenomicRanges_1.54.1       
# [11] GenomeInfoDb_1.38.1         IRanges_2.36.0             
# [13] S4Vectors_0.40.2            BiocGenerics_0.48.1        
# [15] MatrixGenerics_1.14.0       matrixStats_1.2.0          
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.0          viridisLite_0.4.2        
# [3] dplyr_1.1.4               vipor_0.4.5              
# [5] viridis_0.6.4             bitops_1.0-7             
# [7] RCurl_1.98-1.13           RANN_2.6.1               
# [9] bluster_1.11.4            rsvd_1.0.5               
# [11] lifecycle_1.0.4           cluster_2.1.6            
# [13] statmod_1.5.0             magrittr_2.0.3           
# [15] compiler_4.3.2            rlang_1.1.3              
# [17] tools_4.3.2               igraph_2.0.3             
# [19] utf8_1.2.4                S4Arrays_1.2.0           
# [21] dqrng_0.3.2               DelayedArray_0.28.0      
# [23] RColorBrewer_1.1-3        abind_1.4-5              
# [25] BiocParallel_1.36.0       withr_2.5.2              
# [27] grid_4.3.2                fansi_1.0.6              
# [29] beachmat_2.18.0           colorspace_2.1-0         
# [31] edgeR_4.0.3               scales_1.3.0             
# [33] cli_3.6.2                 crayon_1.5.2             
# [35] ragg_1.2.7                generics_0.1.3           
# [37] metapod_1.10.0            DelayedMatrixStats_1.24.0
# [39] ggbeeswarm_0.7.2          zlibbioc_1.48.0          
# [41] parallel_4.3.2            XVector_0.42.0           
# [43] vctrs_0.6.5               Matrix_1.6-4             
# [45] BiocSingular_1.18.0       BiocNeighbors_1.20.2     
# [47] ggrepel_0.9.4             irlba_2.3.5.1            
# [49] beeswarm_0.4.0            systemfonts_1.2.1        
# [51] locfit_1.5-9.8            limma_3.58.1             
# [53] glue_1.7.0                codetools_0.2-19         
# [55] gtable_0.3.4              ScaledMatrix_1.10.0      
# [57] munsell_0.5.0             tibble_3.2.1             
# [59] pillar_1.9.0              GenomeInfoDbData_1.2.11  
# [61] R6_2.5.1                  textshaping_0.3.7        
# [63] sparseMatrixStats_1.14.0  rprojroot_2.0.4          
# [65] lattice_0.22-5            pheatmap_1.0.12          
# [67] Rcpp_1.0.12               gridExtra_2.3            
# [69] SparseArray_1.2.2         pkgconfig_2.0.3  
