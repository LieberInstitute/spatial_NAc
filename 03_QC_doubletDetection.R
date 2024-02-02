#Goal: droplet QC + nuclei QC + doublet detection
#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc

library(SingleCellExperiment)
library(DropletUtils)
library(ggplot2)
library(purrr)
library(here)

#Compile droplet qc information
droplet_paths <- list.files(path = here("processed-data","12_snRNA","droplet_scores"),
                            full.names = TRUE)

names(droplet_paths) <- gsub(x = basename(droplet_paths),
                             pattern = "_droplet_scores.Rdata",
                             replacement = "")

#Read in the droplet scores
e.out <- lapply(droplet_paths, function(x) get(load(x)))

#To make sure we aren't throwing out any cells check if Limited=TRUE and SIG==FALSE
#If both are true, then we could be throwing out non-empty droplets. 
lapply(e.out,function(x){
    table(x$Limited == TRUE & x$FDR>0.001)
})
# $`10c_NAc_SVB`
# 
# FALSE 
# 7178 
# 
# $`11c_NAc_SVB`
# 
# FALSE 
# 4810 
# 
# $`12c_NAc_SVB`
# 
# FALSE 
# 9683 
# 
# $`13c_NAc_SVB`
# 
# FALSE 
# 7729 
# 
# $`14c_NAc_SVB`
# 
# FALSE 
# 8176 
# 
# $`15c_Nac_SVB`
# 
# FALSE 
# 5300 
# 
# $`16c_Nac_SVB`
# 
# FALSE 
# 13984 
# 
# $`17c_Nac_SVB`
# 
# FALSE 
# 7891 
# 
# $`18c_Nac_SVB`
# 
# FALSE 
# 11478 
# 
# $`19c_Nac_SVB`
# 
# FALSE 
# 9640 
# 
# $`1c_NAc_SVB`
# 
# FALSE 
# 8618 
# 
# $`20c_Nac_SVB`
# 
# FALSE 
# 11108 
# 
# $`2c_NAc_SVB`
# 
# FALSE 
# 13779 
# 
# $`3c_NAc_SVB`
# 
# FALSE 
# 8359 
# 
# $`4c_NAc_SVB`
# 
# FALSE 
# 9441 
# 
# $`5c_NAc_SVB`
# 
# FALSE 
# 7851 
# 
# $`6c_NAc_SVB`
# 
# FALSE 
# 16644 
# 
# $`7c_NAc_SVB`
# 
# FALSE 
# 14144 
# 
# $`8c_NAc_SVB`
# 
# FALSE 
# 31312 
# 
# $`9c_NAc_SVB`
# 
# FALSE 
# 9319 

#Also, 
#Another way to look at this
map(e.out, ~ addmargins(table(Signif = .x$FDR <= 0.001, Limited = .x$Limited)))

# $`10c_NAc_SVB`
# Limited
# Signif  FALSE TRUE  Sum
# FALSE   692    0  692
# TRUE    157 6329 6486
# Sum     849 6329 7178
# 
# $`11c_NAc_SVB`
# Limited
# Signif  FALSE TRUE  Sum
# FALSE  1651    0 1651
# TRUE     19 3140 3159
# Sum    1670 3140 4810
# 
# $`12c_NAc_SVB`
# Limited
# Signif  FALSE TRUE  Sum
# FALSE  4557    0 4557
# TRUE     92 5034 5126
# Sum    4649 5034 9683
# 
# $`13c_NAc_SVB`
# Limited
# Signif  FALSE TRUE  Sum
# FALSE  2901    0 2901
# TRUE     48 4780 4828
# Sum    2949 4780 7729
# 
# $`14c_NAc_SVB`
# Limited
# Signif  FALSE TRUE  Sum
# FALSE  1285    0 1285
# TRUE    128 6763 6891
# Sum    1413 6763 8176
# 
# $`15c_Nac_SVB`
# Limited
# Signif  FALSE TRUE  Sum
# FALSE   784    0  784
# TRUE     30 4486 4516
# Sum     814 4486 5300
# 
# $`16c_Nac_SVB`
# Limited
# Signif  FALSE  TRUE   Sum
# FALSE  6544     0  6544
# TRUE    286  7154  7440
# Sum    6830  7154 13984
# 
# $`17c_Nac_SVB`
# Limited
# Signif  FALSE TRUE  Sum
# FALSE  3232    0 3232
# TRUE     69 4590 4659
# Sum    3301 4590 7891
# 
# $`18c_Nac_SVB`
# Limited
# Signif  FALSE  TRUE   Sum
# FALSE  1823     0  1823
# TRUE    177  9478  9655
# Sum    2000  9478 11478
# 
# $`19c_Nac_SVB`
# Limited
# Signif  FALSE TRUE  Sum
# FALSE  2119    0 2119
# TRUE    142 7379 7521
# Sum    2261 7379 9640
# 
# $`1c_NAc_SVB`
# Limited
# Signif  FALSE TRUE  Sum
# FALSE  2460    0 2460
# TRUE    241 5917 6158
# Sum    2701 5917 8618
# 
# $`20c_Nac_SVB`
# Limited
# Signif  FALSE  TRUE   Sum
# FALSE  3958     0  3958
# TRUE    218  6932  7150
# Sum    4176  6932 11108
# 
# $`2c_NAc_SVB`
# Limited
# Signif  FALSE  TRUE   Sum
# FALSE  4799     0  4799
# TRUE    132  8848  8980
# Sum    4931  8848 13779
# 
# $`3c_NAc_SVB`
# Limited
# Signif  FALSE TRUE  Sum
# FALSE  2263    0 2263
# TRUE    337 5759 6096
# Sum    2600 5759 8359
# 
# $`4c_NAc_SVB`
# Limited
# Signif  FALSE TRUE  Sum
# FALSE  2411    0 2411
# TRUE     70 6960 7030
# Sum    2481 6960 9441
# 
# $`5c_NAc_SVB`
# Limited
# Signif  FALSE TRUE  Sum
# FALSE  2478    0 2478
# TRUE     56 5317 5373
# Sum    2534 5317 7851
# 
# $`6c_NAc_SVB`
# Limited
# Signif  FALSE  TRUE   Sum
# FALSE 12646     0 12646
# TRUE     36  3962  3998
# Sum   12682  3962 16644
# 
# $`7c_NAc_SVB`
# Limited
# Signif  FALSE  TRUE   Sum
# FALSE  8104     0  8104
# TRUE    110  5930  6040
# Sum    8214  5930 14144
# 
# $`8c_NAc_SVB`
# Limited
# Signif  FALSE  TRUE   Sum
# FALSE 24005     0 24005
# TRUE    211  7096  7307
# Sum   24216  7096 31312
# 
# $`9c_NAc_SVB`
# Limited
# Signif  FALSE TRUE  Sum
# FALSE  5661    0 5661
# TRUE    181 3477 3658
# Sum    5842 3477 9319


#Not losing any droplets due to number of iterations. Can move on to QC plots and metrics. 






