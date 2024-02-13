#Goal: droplet QC + nuclei QC + doublet detection
#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc
#module load r_nac

library(SingleCellExperiment)
library(DropletUtils)
library(ggplot2)
library(scuttle)
library(scater)
library(purrr)
library(dplyr)
library(tidyr)
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

#Not losing any droplets due to number of iterations.

#Pull knee lower values
#Knee_lowest would work for most but there are a few where we would lose some cells. 
std_out <- readLines(here("code","12_build_sce","02_calculate_droplet_scores.Rout"))
knee_lowers <- std_out[grep("knee_lower =",std_out)][2:21]
knee_lowers <- as.numeric(lapply(strsplit(knee_lowers,split = "="),"[",2))
names(knee_lowers) <- names(e.out)
knee_lowers
# 10c_NAc_SVB 11c_NAc_SVB 12c_NAc_SVB 13c_NAc_SVB 14c_NAc_SVB 15c_Nac_SVB 
# 254         252         214         215         225         234 
# 16c_Nac_SVB 17c_Nac_SVB 18c_Nac_SVB 19c_Nac_SVB  1c_NAc_SVB 20c_Nac_SVB 
# 221         243         294         287         260         256 
# 2c_NAc_SVB  3c_NAc_SVB  4c_NAc_SVB  5c_NAc_SVB  6c_NAc_SVB  7c_NAc_SVB 
# 256         246         266         236         246         271 
# 8c_NAc_SVB  9c_NAc_SVB 
# 265         334

#Create droplet summary table
droplet_summary <- stack(map_int(e.out,nrow)) %>% 
    rename(total_drops=values) %>% 
    left_join(stack(map_int(e.out, ~ sum(.x$FDR < 0.001, na.rm = TRUE)))) %>%
    rename(non_empty=values) %>%
    left_join(stack(knee_lowers)) %>%
    rename(Sample=ind) %>%
    select(Sample,total_drops,non_empty,knee_lower=values)

#Refactor the Sample column so that it is 1-20
droplet_summary$Sample <- factor(x = droplet_summary$Sample,
                                 levels = c("1c_NAc_SVB","2c_NAc_SVB",
                                            "3c_NAc_SVB","4c_NAc_SVB",
                                            "5c_NAc_SVB","6c_NAc_SVB",
                                            "7c_NAc_SVB","8c_NAc_SVB",
                                            "9c_NAc_SVB","10c_NAc_SVB",
                                            "11c_NAc_SVB","12c_NAc_SVB",
                                            "13c_NAc_SVB","14c_NAc_SVB",
                                            "15c_Nac_SVB","16c_Nac_SVB",
                                            "17c_Nac_SVB","18c_Nac_SVB",
                                            "19c_Nac_SVB","20c_Nac_SVB"))

droplet_summary
# Sample total_drops non_empty knee_lower
# 1  10c_NAc_SVB     1516737      6486        254
# 2  11c_NAc_SVB     1126739      3159        252
# 3  12c_NAc_SVB     1674120      5126        214
# 4  13c_NAc_SVB     1301086      4828        215
# 5  14c_NAc_SVB     1432070      6891        225
# 6  15c_Nac_SVB     1322506      4516        234
# 7  16c_Nac_SVB     1680971      7440        221
# 8  17c_Nac_SVB     1313466      4659        243
# 9  18c_Nac_SVB     1623409      9655        294
# 10 19c_Nac_SVB     1375605      7521        287
# 11  1c_NAc_SVB     1284650      6158        260
# 12 20c_Nac_SVB     1631708      7150        256
# 13  2c_NAc_SVB     1903528      8980        256
# 14  3c_NAc_SVB     1368854      6096        246
# 15  4c_NAc_SVB     1568789      7030        266
# 16  5c_NAc_SVB     1219215      5373        236
# 17  6c_NAc_SVB      984751      3998        246
# 18  7c_NAc_SVB     1371545      6040        271
# 19  8c_NAc_SVB     1693619      7307        265
# 20  9c_NAc_SVB      876588      3658        334

#Write out file. 
write.csv(x = droplet_summary,
          file = here("processed-data","12_snRNA","NAc_snRNA_droplet_summary.csv"),
          row.names = FALSE,
          quote = FALSE)

#Make a barplot summarizing the number of empty and non-empty droplets. 
droplet_barplot <- droplet_summary %>%
    mutate(empty = total_drops - non_empty) %>%
    select(-total_drops) %>%
    select(-knee_lower) %>%
    pivot_longer(!Sample,names_to = "drop_type",values_to = "number_drops") %>%
    ggplot(aes(x = Sample,y=number_drops,fill = drop_type)) +
    geom_col() +
    scale_y_continuous(trans = "log10") +
    labs(x = "Sample",
         y = "Number of Droplets",
         fill = "Droplet Status") +
    theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggsave(plot = droplet_barplot,filename = here("plots","12_snRNA","droplet_barplot_per_sample.png"))


#Load in the sce object
load(here("processed-data","12_snRNA","sce_raw.rds"),verbose = TRUE)
# Loading objects:
#     sce

dim(sce)
#[1]    36601 28269956

#Refactor the sample names so QC plot x-axes are in correct order. 
sce$Sample <- factor(x = sce$Sample,
                     levels = c("1c_NAc_SVB","2c_NAc_SVB",
                                "3c_NAc_SVB","4c_NAc_SVB",
                                "5c_NAc_SVB","6c_NAc_SVB",
                                "7c_NAc_SVB","8c_NAc_SVB",
                                "9c_NAc_SVB","10c_NAc_SVB",
                                "11c_NAc_SVB","12c_NAc_SVB",
                                "13c_NAc_SVB","14c_NAc_SVB",
                                "15c_Nac_SVB","16c_Nac_SVB",
                                "17c_Nac_SVB","18c_Nac_SVB",
                                "19c_Nac_SVB","20c_Nac_SVB"))

#### Eliminate empty droplets ####
e.out.all <- do.call("rbind", e.out)[colnames(sce), ]
sce <- sce[, which(e.out.all$FDR <= 0.001)]

dim(sce)
#[1]  36601 122071

#Save object
save(sce,file = here("processed-data","12_snRNA","sce_emptyDrops_removed.rds"))

####Begin QC
sce <- scuttle::addPerCellQC(sce,subsets = list(Mito=which(seqnames(sce) == "chrM")))

#Plot mitochondria vs detected
#All samples together
mito_vs_detected <- plotColData(object = sce,
                                y = "subsets_Mito_percent",
                                x = "detected",
                                colour_by = "Sample")

ggsave(filename = here("plots","12_snRNA","nuclei_QC","mito","mito_vs_detected_bySample.png"),
       plot = mito_vs_detected)

#Now samples separately.
for(i in unique(sce$Sample)){
    print(i)
    x <- plotColData(object = sce[,sce$Sample == i],
                     y = "subsets_Mito_percent",
                     x = "detected",
                     colour_by = "Sample")
    ggsave(filename = here("plots","12_snRNA","nuclei_QC","mito",
                           paste0("mito_vs_detected_",i,"_only.png")),
           plot = x)
}

#########################################
############### High mito ###############
#########################################
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher", batch = sce$Sample)
table(sce$Sample,sce$high_mito)
# FALSE TRUE
# 1c_NAc_SVB   5555  603
# 2c_NAc_SVB   8265  715
# 3c_NAc_SVB   5462  634
# 4c_NAc_SVB   6489  541
# 5c_NAc_SVB   4875  498
# 6c_NAc_SVB   3630  368
# 7c_NAc_SVB   5140  900
# 8c_NAc_SVB   7055  252
# 9c_NAc_SVB   3506  152
# 10c_NAc_SVB  6063  423
# 11c_NAc_SVB  2951  208
# 12c_NAc_SVB  4637  489
# 13c_NAc_SVB  4336  492
# 14c_NAc_SVB  6361  530
# 15c_Nac_SVB  4132  384
# 16c_Nac_SVB  7027  413
# 17c_Nac_SVB  4189  470
# 18c_Nac_SVB  8960  695
# 19c_Nac_SVB  6937  584
# 20c_Nac_SVB  6313  837

nMAD3_vln <- plotColData(sce,x = "Sample",y= "subsets_Mito_percent",colour_by = "high_mito") +
    ggtitle(paste0("Mito Percent\n3 MADs")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 5,lty = 2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = nMAD3_vln,filename = here("plots","12_snRNA",
                                        "nuclei_QC","mito",
                                        "nMAD3_Mito_Pct_violin.png"))

#Simialr to the Comment in https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/51d15ef9f5f2c4c53f55e22e3fe467de1a724668/10x_all-FACS-n10_2021rev_step01_processing-QC_MNT.R#L4
#MAD approach is unnecessarily throwing out cells because the distribution is centered around 0

#Check the number of nuclei 
sce$high_mito_pct <- ifelse(sce$subsets_Mito_percent > 5.0,
                            TRUE,
                            FALSE)

table(sce$high_mito_pct)
# FALSE   TRUE 
# 121825    246 

mito_pct_vln <- plotColData(sce,x = "Sample",y= "subsets_Mito_percent",colour_by = "high_mito_pct") +
    ggtitle(paste0("Mito Percent\n5% cutoff")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 5,lty = 2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = mito_pct_vln,filename = here("plots","12_snRNA","nuclei_QC","mito","FivePct_Mito_Pct_violin.png"))


#########################################
########## Low Library Size #############
#########################################
# ## low library size
#Plot library size per sample to look at distributions. 
lib_size_violin <- plotColData(sce, x = "Sample", y = "sum",colour_by = "Sample") +
    scale_y_log10() +
    ggtitle("Total UMIs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = lib_size_violin,filename = here("plots","12_snRNA","nuclei_QC","library_size","library_size_violin.png"))

#Some with unimodal and some with bimodal distributions. 
#Try nMADs=3
##3
sce$low_lib_3 <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$Sample,nmads = 3)
table(sce$Sample,sce$low_lib_3)


nMAD3_lib_vln <- plotColData(sce, x = "Sample", y = "sum",colour_by = "low_lib_3") +
    scale_y_log10() +
    ggtitle("Total UMIs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = nMAD3_lib_vln,filename = here("plots","12_snRNA","nuclei_QC","nMAD3_library_size_violin.png"))

#Some are bimodal. Some are unimodal. Doesn't track with sort type. 
#PI sort should have both neurons and glia, while PI_NeuN should have primarily neurons. 
#For some sort days it seems to track, for others it doesn't. Is it due to the staining success from each data?
#snRNA_data should actually be snRNA_data. 
lib_size_violin <- plotColData(sce, x = "Sample", y = "sum",colour_by = "snRNA_data") +
    scale_y_log10() +
    ggtitle("Total UMIs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = lib_size_violin,filename = here("plots","12_snRNA","nuclei_QC","library_size_violin_coloredbydate.png"))

#While a totally different metric, we expect that neurons and glia should have different numbers of genes/nucleus.
#Plotting the genes might make it easier to see any sort of pattern. 
detected_violin <- plotColData(sce, x = "Sample", y = "detected",colour_by = "Sort") +
    scale_y_log10() +
    ggtitle("Total UMIs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = detected_violin,filename = here("plots","12_snRNA","nuclei_QC","number_of_genes_violin_coloredbySort.png"))

#Date as well 
detected_violin <- plotColData(sce, x = "Sample", y = "detected",colour_by = "snRNA_data") +
    scale_y_log10() +
    ggtitle("Total UMIs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = detected_violin,filename = here("plots","12_snRNA","nuclei_QC","number_of_genes_violin_coloredbySortDate.png"))


#calculate nMAD by sort date. 
#Reasoning here is that staining success could change by batch date. 
#First change the sort name. 
colnames(colData(sce))[6] <- "snRNA_date"

# sce$low_lib_date_2 <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$snRNA_date,nmads = 2)
# sce$low_lib_date_3 <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$snRNA_date,nmads = 3)
# 
# plotColData(sce, x = "Sample", y = "detected",colour_by = "low_lib_date_2") +
#     scale_y_log10() +
#     ggtitle("Total UMIs") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# plotColData(sce, x = "Sample", y = "detected",colour_by = "low_lib_date_3") +
#     scale_y_log10() +
#     ggtitle("Total UMIs") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# sce$low_lib_sort <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$Sort,nmads = 2)
# plotColData(sce, x = "Sample", y = "detected",colour_by = "low_lib_sort") +
#     scale_y_log10() +
#     ggtitle("Total UMIs") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))

###make some basic QC cutoffs and move on to see if any PI_NeuN samples have a substantial amount of glia 
#Will do total umis < 600,detected genes > 300, and mitochondrial percentage > 50
qc_lib_600 <- sce$sum < 600
qc_genes_500 <- sce$detected < 500
qc_mito_5 <- sce$subsets_Mito_percent > 5

sce$discard_basic <- qc_lib_600 | qc_genes_500 | qc_mito_5
# FALSE   TRUE 
# 118358   3713 

qc_t <- addmargins(table(sce$Sample, sce$discard_basic))

qc_t








