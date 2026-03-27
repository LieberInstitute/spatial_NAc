# cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
# module load r_nac
library(SpatialExperiment)
library(here)

spe <- readRDS(here("processed-data","xenium","SPEs","spe_NormCounts.Rds"))

dim(spe)

xenium_genes <- rowData(spe)[rowData(spe)$Type == "Gene Expression",]

xenium_df <- as.data.frame(xenium_genes[,c("Symbol","ID")])
colnames(xenium_df)[2] <- "Ensembl_ID"
rownames(xenium_df) <- NULL


write.csv(x = xenium_df,file = here("processed-data","12_snRNA",
                                    "xenium_panel_supptable17.csv"))


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
#   [1] here_1.0.1                  SpatialExperiment_1.12.0   
# [3] SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
# [5] Biobase_2.62.0              GenomicRanges_1.54.1       
# [7] GenomeInfoDb_1.38.1         IRanges_2.36.0             
# [9] S4Vectors_0.40.2            BiocGenerics_0.48.1        
# [11] MatrixGenerics_1.14.0       matrixStats_1.2.0          
# 
# loaded via a namespace (and not attached):
#   [1] crayon_1.5.2            magick_2.8.1            DelayedArray_0.28.0    
# [4] rjson_0.2.21            RCurl_1.98-1.13         rprojroot_2.0.4        
# [7] grid_4.3.2              abind_1.4-5             bitops_1.0-7           
# [10] compiler_4.3.2          Rcpp_1.0.12             XVector_0.42.0         
# [13] lattice_0.22-5          SparseArray_1.2.2       GenomeInfoDbData_1.2.11
# [16] magrittr_2.0.3          Matrix_1.6-4            tools_4.3.2            
# [19] zlibbioc_1.48.0         S4Arrays_1.2.0         
