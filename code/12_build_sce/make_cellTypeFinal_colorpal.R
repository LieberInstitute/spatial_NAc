#Completed in interactive job

library(Polychrome)
library(here)

cluster_cols <- Polychrome::createPalette(N = 21,
                                          seedcolors = c("#000000", "#E69F00",
                                                         "#56B4E9", "#009E73", 
                                                         "#F0E442", "#0072B2", 
                                                         "#D55E00", "#CC79A7"),
                                          target = "normal")

cluster_cols
NC1       NC2       NC3       NC4       NC5       NC6       NC7       NC8 
"#4F4753" "#ECA31C" "#58B6ED" "#0D9F72" "#F2E642" "#0077B9" "#D95F00" "#D079AA" 
NC9      NC10      NC11      NC12      NC13      NC14      NC15      NC16 
"#D00DFF" "#35FB00" "#F80091" "#FF0016" "#2A4BF9" "#F0DEC4" "#FB3DD9" "#7A0096" 
NC17      NC18      NC19      NC20      NC21 
"#854222" "#A7F281" "#0DFBFA" "#5C6300" "#E2D4FB" 

#colors 14 and 21 are very light. 
cluster_cols[21] <- "black"
cluster_cols[14] <- "grey"

save(cluster_cols,file = here("processed-data","12_snRNA","070724_21colors_celltypeFinal.rda"))

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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] here_1.0.1       Polychrome_1.5.1
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_2.1-0     compiler_4.3.2       rprojroot_2.0.4     
# [4] scatterplot3d_0.3-44
