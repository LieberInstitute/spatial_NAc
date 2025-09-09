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


#GJA1 featureplot
gja1_Feature <- plotReducedDim(object = sce,
                             dimred = "tSNE_HARMONY",
                             color_by = "GJA1",
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
ggsave(plot = gja1_Feature,
       filename = here("plots","12_snRNA","Supplementary","Additional_Feature_Plots","GJA1_FeaturePlot_legend.pdf"),
       height = 8,
       width = 8)
#without legend
ggsave(plot = gja1_Feature + theme(legend.position = "none"),
       filename = here("plots","12_snRNA","Supplementary","Additional_Feature_Plots","GJA1_FeaturePlot_nolegend.png"),
       height = 8,
       width = 8)


#MOBP featureplot
MOBP_Feature <- plotReducedDim(object = sce,
                               dimred = "tSNE_HARMONY",
                               color_by = "MOBP",
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
ggsave(plot = MOBP_Feature,
       filename = here("plots","12_snRNA","Supplementary","Additional_Feature_Plots","MOBP_FeaturePlot_legend.pdf"),
       height = 8,
       width = 8)
#without legend
ggsave(plot = MOBP_Feature + theme(legend.position = "none"),
       filename = here("plots","12_snRNA","Supplementary","Additional_Feature_Plots","MOBP_FeaturePlot_nolegend.png"),
       height = 8,
       width = 8)


#PDGFRA featureplot
PDGFRA_Feature <- plotReducedDim(object = sce,
                               dimred = "tSNE_HARMONY",
                               color_by = "PDGFRA",
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
ggsave(plot = PDGFRA_Feature,
       filename = here("plots","12_snRNA","Supplementary","Additional_Feature_Plots","PDGFRA_FeaturePlot_legend.pdf"),
       height = 8,
       width = 8)
#without legend
ggsave(plot = PDGFRA_Feature + theme(legend.position = "none"),
       filename = here("plots","12_snRNA","Supplementary","Additional_Feature_Plots","PDGFRA_FeaturePlot_nolegend.png"),
       height = 8,
       width = 8)



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
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     
# 
# other attached packages:
#   [1] here_1.0.1                  scran_1.30.0               
# [3] dplyr_1.1.4                 scater_1.30.1              
# [5] scuttle_1.12.0              ggplot2_3.5.1              
# [7] sessioninfo_1.2.2           ComplexHeatmap_2.18.0      
# [9] SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
# [11] Biobase_2.62.0              GenomicRanges_1.54.1       
# [13] GenomeInfoDb_1.38.1         IRanges_2.36.0             
# [15] S4Vectors_0.40.2            BiocGenerics_0.48.1        
# [17] MatrixGenerics_1.14.0       matrixStats_1.2.0          
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.0          viridisLite_0.4.2        
# [3] farver_2.1.1              vipor_0.4.5              
# [5] viridis_0.6.4             bitops_1.0-7             
# [7] RCurl_1.98-1.13           bluster_1.11.4           
# [9] digest_0.6.33             rsvd_1.0.5               
# [11] lifecycle_1.0.4           cluster_2.1.6            
# [13] statmod_1.5.0             magrittr_2.0.3           
# [15] compiler_4.3.2            rlang_1.1.3              
# [17] tools_4.3.2               igraph_2.0.3             
# [19] utf8_1.2.4                labeling_0.4.3           
# [21] S4Arrays_1.2.0            dqrng_0.3.2              
# [23] DelayedArray_0.28.0       RColorBrewer_1.1-3       
# [25] abind_1.4-5               BiocParallel_1.36.0      
# [27] withr_2.5.2               fansi_1.0.6              
# [29] beachmat_2.18.0           colorspace_2.1-0         
# [31] edgeR_4.0.3               scales_1.3.0             
# [33] iterators_1.0.14          cli_3.6.2                
# [35] crayon_1.5.2              ragg_1.2.7               
# [37] generics_0.1.3            metapod_1.10.0           
# [39] rjson_0.2.21              DelayedMatrixStats_1.24.0
# [41] ggbeeswarm_0.7.2          zlibbioc_1.48.0          
# [43] parallel_4.3.2            XVector_0.42.0           
# [45] vctrs_0.6.5               Matrix_1.6-4             
# [47] BiocSingular_1.18.0       GetoptLong_1.0.5         
# [49] BiocNeighbors_1.20.2      ggrepel_0.9.4            
# [51] irlba_2.3.5.1             clue_0.3-65              
# [53] beeswarm_0.4.0            systemfonts_1.2.1        
# [55] locfit_1.5-9.8            foreach_1.5.2            
# [57] limma_3.58.1              glue_1.7.0               
# [59] codetools_0.2-19          cowplot_1.1.1            
# [61] shape_1.4.6               gtable_0.3.4             
# [63] ScaledMatrix_1.10.0       munsell_0.5.0            
# [65] tibble_3.2.1              pillar_1.9.0             
# [67] GenomeInfoDbData_1.2.11   circlize_0.4.16          
# [69] R6_2.5.1                  textshaping_0.3.7        
# [71] sparseMatrixStats_1.14.0  rprojroot_2.0.4          
# [73] doParallel_1.0.17         lattice_0.22-5           
# [75] png_0.1-8                 Rcpp_1.0.12              
# [77] gridExtra_2.3             SparseArray_1.2.2        
# [79] pkgconfig_2.0.3           GlobalOptions_0.1.2      
# 
# 
# 
# 
# 
