#Goal: droplet QC + nuclei QC + doublet detection
#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc
#module load r_nac

library(SingleCellExperiment)
library(DropletUtils)
library(scDblFinder)
library(HDF5Array)
library(rafalib)
library(ggplot2)
library(scuttle)
library(scran)
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
#Batch with snRNA_date
#name is currently snRNA_data --> Need to change to snRNA_date
colnames(colData(sce))[6] <- "snRNA_date"
#Try nMADs=3
##3
sce$low_lib_3 <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$snRNA_date,nmads = 3)
table(sce$Sample,sce$low_lib_3)
#             FALSE TRUE
# 1c_NAc_SVB   6024  134
# 2c_NAc_SVB   8830  150
# 3c_NAc_SVB   6070   26
# 4c_NAc_SVB   6995   35
# 5c_NAc_SVB   5306   67
# 6c_NAc_SVB   3971   27
# 7c_NAc_SVB   5904  136
# 8c_NAc_SVB   7017  290
# 9c_NAc_SVB   3658    0
# 10c_NAc_SVB  6474   12
# 11c_NAc_SVB  3151    8
# 12c_NAc_SVB  5093   33
# 13c_NAc_SVB  4800   28
# 14c_NAc_SVB  6855   36
# 15c_Nac_SVB  4516    0
# 16c_Nac_SVB  7440    0
# 17c_Nac_SVB  4659    0
# 18c_Nac_SVB  9655    0
# 19c_Nac_SVB  7521    0
# 20c_Nac_SVB  7150    0

nMAD3_lib_vln <- plotColData(sce, x = "Sample", y = "sum",colour_by = "low_lib_3") +
    scale_y_log10() +
    ggtitle("Total UMIs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = nMAD3_lib_vln,filename = here("plots","12_snRNA","nuclei_QC","library_size","nMAD3_library_size_violin.png"))

#nMAD3 is not removing enough low quality cells. Cells with <1000 UMIs being kept and no cells from last batch being kept. 
sce$low_lib_2 <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$snRNA_date,nmads = 2)
table(sce$Sample,sce$low_lib_2)
#             FALSE TRUE
# 1c_NAc_SVB   5604  554
# 2c_NAc_SVB   8452  528
# 3c_NAc_SVB   5929  167
# 4c_NAc_SVB   6875  155
# 5c_NAc_SVB   5033  340
# 6c_NAc_SVB   3846  152
# 7c_NAc_SVB   5721  319
# 8c_NAc_SVB   6608  699
# 9c_NAc_SVB   3326  332
# 10c_NAc_SVB  6338  148
# 11c_NAc_SVB  3112   47
# 12c_NAc_SVB  4803  323
# 13c_NAc_SVB  4625  203
# 14c_NAc_SVB  6572  319
# 15c_Nac_SVB  4399  117
# 16c_Nac_SVB  7252  188
# 17c_Nac_SVB  4575   84
# 18c_Nac_SVB  9173  482
# 19c_Nac_SVB  7199  322
# 20c_Nac_SVB  6866  284

#Just how many cells are removed
table(sce$low_lib_2)
# FALSE   TRUE 
# 116308   5763

nMAD2_lib_vln <- plotColData(sce, x = "Sample", y = "sum",colour_by = "low_lib_2") +
    scale_y_log10() +
    ggtitle("Total UMIs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = nMAD3_lib_vln,filename = here("plots","12_snRNA","nuclei_QC","library_size","nMAD2_library_size_violin.png"))


#Test for unimodality with Hartigans’ Dip Test for Unimodality
#Null hypothesis is that distribution is unimodal
#empty dataframe to input results into. 
sum_res_df <- data.frame(Sample = unique(sce$Sample),
                         p = NA,
                         row.names = unique(sce$Sample))

for(i in unique(sce$Sample)){
    x <- subset(colData(sce),subset=(Sample == i))$sum
    dip.test.results <- dip.test(x)
    sum_res_df[i,"p"] <- dip.test.results$p
}

#Adjust the p-value
sum_res_df$is.unimodal <- ifelse(sum_res_df$p <= 0.05,
                                 FALSE,
                                 TRUE)

colData(sce) <- merge(x = colData(sce),
                      y = sum_res_df,
                      by = "Sample")

plotColData(sce, x = "Sample", y = "sum",colour_by = "is.unimodal") +
    scale_y_log10() +
    ggtitle("Total UMIs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Some are bimodal. Some are unimodal. Doesn't track with sort type. 
#PI sort should have both neurons and glia, while PI_NeuN should have primarily neurons. 
#For some sort days it seems to track, for others it doesn't. Is it due to the staining success from each data?
#snRNA_data should actually be snRNA_data. 
#First change the sort name. 
colnames(colData(sce))[6] <- "snRNA_date"
lib_size_violin <- plotColData(sce, x = "Sample", y = "sum",colour_by = "snRNA_date") +
    scale_y_log10() +
    ggtitle("Total UMIs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = lib_size_violin,filename = here("plots","12_snRNA","nuclei_QC","library_size","library_size_violin_coloredbydate.png"))

#While a totally different metric, we expect that neurons and glia should have different numbers of genes/nucleus.
#Plotting the genes might make it easier to see any sort of pattern. 
detected_violin <- plotColData(sce, x = "Sample", y = "detected",colour_by = "Sort") +
    scale_y_log10() +
    ggtitle("Total UMIs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = detected_violin,filename = here("plots","12_snRNA","nuclei_QC","number_genes","number_of_genes_violin_coloredbySort.png"))

#Date as well 
detected_violin <- plotColData(sce, x = "Sample", y = "detected",colour_by = "snRNA_date") +
    scale_y_log10() +
    ggtitle("Total UMIs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = detected_violin,filename = here("plots","12_snRNA","nuclei_QC","number_genes","number_of_genes_violin_coloredbySortDate.png"))


###make some basic QC cutoffs and move on to see if any PI_NeuN samples have a substantial amount of glia 
#Will do total umis < 600,detected genes > 300, and mitochondrial percentage > 50
qc_lib_600 <- sce$sum < 600
qc_genes_500 <- sce$detected < 500
qc_mito_5 <- sce$subsets_Mito_percent > 5

sce$discard_basic <- qc_lib_600 | qc_genes_500 | qc_mito_5


qc_t <- addmargins(table(sce$Sample, sce$discard_basic))

qc_t
#              FALSE   TRUE    Sum
# 1c_NAc_SVB    5868    290   6158
# 2c_NAc_SVB    8732    248   8980
# 3c_NAc_SVB    6036     60   6096
# 4c_NAc_SVB    6970     60   7030
# 5c_NAc_SVB    5276     97   5373
# 6c_NAc_SVB    3950     48   3998
# 7c_NAc_SVB    5838    202   6040
# 8c_NAc_SVB    6859    448   7307
# 9c_NAc_SVB    3417    241   3658
# 10c_NAc_SVB   6413     73   6486
# 11c_NAc_SVB   3129     30   3159
# 12c_NAc_SVB   4953    173   5126
# 13c_NAc_SVB   4723    105   4828
# 14c_NAc_SVB   6721    170   6891
# 15c_Nac_SVB   4398    118   4516
# 16c_Nac_SVB   7255    185   7440
# 17c_Nac_SVB   4573     86   4659
# 18c_Nac_SVB   9190    465   9655
# 19c_Nac_SVB   7202    319   7521
# 20c_Nac_SVB   6855    295   7150
# Sum         118358   3713 122071


#### Doublet detection ####
## To speed up, run on sample-level top-HVGs - just take top 1000
set.seed(1234)

colData(sce)$doubletScore <- NA

for (i in splitit(sce$Sample)) {
    sce_temp <- sce[, i]
    ## To speed up, run on sample-level top-HVGs - just take top 1000
    normd <- logNormCounts(sce_temp)
    geneVar <- modelGeneVar(normd)
    topHVGs <- getTopHVGs(geneVar, n = 1000)
    
    dbl_dens <- computeDoubletDensity(normd, subset.row = topHVGs)
    colData(sce)$doubletScore[i] <- dbl_dens
}

summary(sce$doubletScore)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00000  0.07522  0.23172  0.53328  0.50616 34.79662 

## Visualize doublet scores ##
dbl_df <- colData(sce) %>%
    as.data.frame() %>%
    select(Sample, doubletScore)

dbl_box_plot <- dbl_df %>%
    ggplot(aes(x = Sample, y = doubletScore, fill = Sample)) +
    geom_boxplot() +
    labs(x = "Sample") +
    geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
    coord_flip() +
    theme_bw()

ggsave(dbl_box_plot, filename = here("plots","12_snRNA",
                                     "nuclei_QC","doublet_score",
                                     "doublet_scores_boxplot.png"))

dbl_density_plot <- dbl_df %>%
    ggplot(aes(x = doubletScore,fill = Sample)) +
    geom_density() +
    labs(x = "doublet score") +
    theme_bw()

ggsave(dbl_density_plot, filename = here("plots","12_snRNA",
                                         "nuclei_QC","doublet_score",
                                         "doublet_scores_desnity.png"))

dbl_df %>%
    group_by(Sample) %>%
    summarize(
        median = median(doubletScore),
        q95 = quantile(doubletScore, .95),
        drop = sum(doubletScore >= 5),
        drop_precent = 100 * drop / n()
    )
# # A tibble: 20 × 5
# Sample      median   q95  drop drop_precent
# <fct>        <dbl> <dbl> <int>        <dbl>
#     1 1c_NAc_SVB   0.209 1.21     70        1.14 
# 2 2c_NAc_SVB   0.198 1.20    154        1.71 
# 3 3c_NAc_SVB   0.158 1.39    119        1.95 
# 4 4c_NAc_SVB   0.211 2.28    182        2.59 
# 5 5c_NAc_SVB   0.204 1.10     70        1.30 
# 6 6c_NAc_SVB   0.280 1.25     64        1.60 
# 7 7c_NAc_SVB   0.230 1.24    101        1.67 
# 8 8c_NAc_SVB   0.234 0.935    92        1.26 
# 9 9c_NAc_SVB   0.490 1.51     35        0.957
# 10 10c_NAc_SVB  0.195 1.71    149        2.30 
# 11 11c_NAc_SVB  0.461 1.85     51        1.61 
# 12 12c_NAc_SVB  0.267 1.31     54        1.05 
# 13 13c_NAc_SVB  0.338 1.59     57        1.18 
# 14 14c_NAc_SVB  0.179 1.35    117        1.70 
# 15 15c_Nac_SVB  0.280 1.16     30        0.664
# 16 16c_Nac_SVB  0.164 1.47    143        1.92 
# 17 17c_Nac_SVB  0.298 1.71     78        1.67 
# 18 18c_Nac_SVB  0.212 1.47    184        1.91 
# 19 19c_Nac_SVB  0.211 1.55     79        1.05 
# 20 20c_Nac_SVB  0.229 1.20     85        1.19

table(sce$discard_basic, sce$doubletScore >= 5)
#        FALSE   TRUE
# FALSE 116444   1914
# TRUE    3713      0
#No cells that are discarded by basic QC are also discarded by doubletscore

#Save object
saveHDF5SummarizedExperiment(sce,
                             here("processed-data","12_snRNA",
                                  "sce_emptyDrops_removed_withQC"),
                             replace = TRUE)

#Remove the cells that don't meet the basic QC cutoffs 
sce <- sce[,!sce$discard_basic]
dim(sce)
#[1]  36601 118358

#Save object
saveHDF5SummarizedExperiment(sce,
                             here("processed-data","12_snRNA",
                                  "sce_clean"),
                             replace = TRUE)


print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
# [1] "Reproducibility information:"
# [1] "2024-02-13 19:46:41 EST"
# user   system  elapsed 
# 744.133    4.288 1407.265 
# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-02-13
# pandoc   3.1.3 @ /jhpce/shared/libd/core/r_nac/1.0/nac_env/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                * 1.4-5     2016-07-21 [1] CRAN (R 4.3.2)
# beachmat               2.18.0    2023-10-24 [1] Bioconductor
# beeswarm               0.4.0     2021-06-01 [1] CRAN (R 4.3.2)
# Biobase              * 2.62.0    2023-10-24 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-01 [1] Bioconductor
# BiocIO                 1.12.0    2023-10-24 [1] Bioconductor
# BiocNeighbors          1.20.2    2024-01-07 [1] Bioconductor 3.18 (R 4.3.2)
# BiocParallel           1.36.0    2023-10-24 [1] Bioconductor
# BiocSingular           1.18.0    2023-10-24 [1] Bioconductor
# Biostrings             2.70.1    2023-10-25 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.2)
# bluster                1.11.4    2024-02-02 [1] Github (LTLA/bluster@17dd9c8)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.2)
# cluster                2.1.6     2023-12-01 [1] CRAN (R 4.3.2)
# codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.0)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# data.table             1.14.10   2023-12-08 [1] CRAN (R 4.3.2)
# DelayedArray         * 0.28.0    2023-10-24 [1] Bioconductor
# DelayedMatrixStats     1.24.0    2023-10-24 [1] Bioconductor
# dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.3.2)
# dqrng                  0.3.2     2023-11-29 [1] CRAN (R 4.3.2)
# DropletUtils         * 1.22.0    2023-10-24 [1] Bioconductor
# edgeR                  4.0.3     2023-12-10 [1] Bioconductor 3.18 (R 4.3.2)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.2)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.38.1    2023-11-08 [1] Bioconductor
# GenomeInfoDbData       1.2.11    2023-12-12 [1] Bioconductor
# GenomicAlignments      1.38.0    2023-10-24 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-29 [1] Bioconductor
# ggbeeswarm             0.7.2     2023-04-29 [1] CRAN (R 4.3.2)
# ggplot2              * 3.4.4     2023-10-12 [1] CRAN (R 4.3.1)
# ggrepel                0.9.4     2023-10-13 [1] CRAN (R 4.3.2)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.2)
# gridExtra              2.3       2017-09-09 [1] CRAN (R 4.3.2)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.1)
# HDF5Array            * 1.30.0    2023-10-24 [1] Bioconductor
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.2)
# igraph                 1.6.0     2023-12-11 [1] CRAN (R 4.3.2)
# IRanges              * 2.36.0    2023-10-24 [1] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.3.2)
# jsonlite               1.8.8     2023-12-04 [1] CRAN (R 4.3.2)
# lattice                0.22-5    2023-10-24 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.2)
# limma                  3.58.1    2023-10-31 [1] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [1] CRAN (R 4.3.2)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# MASS                   7.3-60    2023-05-04 [1] CRAN (R 4.3.0)
# Matrix               * 1.6-4     2023-11-30 [1] CRAN (R 4.3.2)
# MatrixGenerics       * 1.14.0    2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.2)
# metapod                1.10.0    2023-10-24 [1] Bioconductor
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# purrr                * 1.0.2     2023-08-10 [1] CRAN (R 4.3.1)
# R.methodsS3            1.8.2     2022-06-13 [1] CRAN (R 4.3.2)
# R.oo                   1.26.0    2024-01-24 [1] CRAN (R 4.3.2)
# R.utils                2.12.3    2023-11-18 [1] CRAN (R 4.3.2)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.3.2)
# RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.2)
# RCurl                  1.98-1.13 2023-11-02 [1] CRAN (R 4.3.2)
# restfulr               0.0.15    2022-06-16 [1] CRAN (R 4.3.2)
# rhdf5                * 2.46.1    2023-11-29 [1] Bioconductor 3.18 (R 4.3.2)
# rhdf5filters           1.14.1    2023-11-06 [1] Bioconductor
# Rhdf5lib               1.24.1    2023-12-11 [1] Bioconductor 3.18 (R 4.3.2)
# rjson                  0.2.21    2022-01-09 [1] CRAN (R 4.3.2)
# rlang                  1.1.3     2024-01-10 [1] CRAN (R 4.3.2)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.2)
# Rsamtools              2.18.0    2023-10-24 [1] Bioconductor
# rsvd                   1.0.5     2021-04-16 [1] CRAN (R 4.3.2)
# rtracklayer            1.62.0    2023-10-24 [1] Bioconductor
# S4Arrays             * 1.2.0     2023-10-24 [1] Bioconductor
# S4Vectors            * 0.40.2    2023-11-23 [1] Bioconductor 3.18 (R 4.3.2)
# ScaledMatrix           1.10.0    2023-10-24 [1] Bioconductor
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.2)
# scater               * 1.30.1    2023-11-16 [1] Bioconductor
# scDblFinder          * 1.16.0    2023-10-24 [1] Bioconductor
# scran                * 1.30.0    2023-10-24 [1] Bioconductor
# scuttle              * 1.12.0    2023-10-24 [1] Bioconductor
# sessioninfo            1.2.2     2021-12-06 [1] CRAN (R 4.3.2)
# SingleCellExperiment * 1.24.0    2023-10-24 [1] Bioconductor
# SparseArray          * 1.2.2     2023-11-07 [1] Bioconductor
# sparseMatrixStats      1.14.0    2023-10-24 [1] Bioconductor
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.2)
# SummarizedExperiment * 1.32.0    2023-10-24 [1] Bioconductor
# tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                * 1.3.0     2023-01-24 [1] CRAN (R 4.3.0)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.2)
# vipor                  0.4.5     2017-03-22 [1] CRAN (R 4.3.2)
# viridis                0.6.4     2023-07-22 [1] CRAN (R 4.3.2)
# viridisLite            0.4.2     2023-05-02 [1] CRAN (R 4.3.0)
# withr                  2.5.2     2023-10-30 [1] CRAN (R 4.3.1)
# xgboost                1.7.6.1   2023-12-06 [1] CRAN (R 4.3.2)
# XML                    3.99-0.16 2023-11-29 [1] CRAN (R 4.3.2)
# XVector                0.42.0    2023-10-24 [1] Bioconductor
# yaml                   2.3.8     2023-12-11 [1] CRAN (R 4.3.2)
# zlibbioc               1.48.0    2023-10-24 [1] Bioconductor
# 
# [1] /jhpce/shared/libd/core/r_nac/1.0/nac_env/lib/R/library
# 
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
