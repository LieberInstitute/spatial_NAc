library(SingleCellExperiment)
library(DropletUtils)
library(scDblFinder)
library(HDF5Array)
library(rafalib)
library(ggplot2)
library(scuttle)
library(readr)
library(scran)
library(scater)
library(purrr)
library(dplyr)
library(tidyr)
library(here)

#Compile droplet qc information
print("Reading in droplet paths")
droplet_paths <- list.files(path = here("processed-data","21_Reiner_snRNAseq","droplet_scores"),
                            full.names = TRUE)

names(droplet_paths) <- gsub(x = basename(droplet_paths),
                             pattern = "_droplet_scores.Rdata",
                             replacement = "")

#Read in the droplet scores
e.out <- lapply(droplet_paths, function(x) get(load(x)))
#To make sure we aren't throwing out any cells check if Limited=TRUE and SIG==FALSE
#If both are true, then we could be throwing out non-empty droplets. 
map(e.out, ~ addmargins(table(Signif = .x$FDR <= 0.001, Limited = .x$Limited)))

#Pull knee lower values
std_out <- list.files(here("processed-data", "21_Reiner_snRNAseq", "logs"), pattern = "emptydrops", full.names = TRUE)
std_out <- map(std_out, readLines)

#118th line of each iteration of the list is the knee_lower value
knee_lower <- as.character(lapply(X = std_out,"[",118))
knee_lower

#Just keep the values
knee_lower <- as.numeric(lapply(strsplit(knee_lower,split = "="),"[",2))

#151st element of the std_out list contains the sample name
sample_name <- as.character(lapply(std_out,"[",109))
sample_name

sample_name <- as.character(lapply(strsplit(sample_name,split = " "),"[",3))

names(knee_lower) <- sample_name

#Create droplet summary table
droplet_summary <- stack(map_int(e.out,nrow)) %>% 
    rename(total_drops=values) %>% 
    left_join(stack(map_int(e.out, ~ sum(.x$FDR < 0.001, na.rm = TRUE)))) %>%
    rename(non_empty=values) %>%
    left_join(stack(knee_lower)) %>%
    rename(Sample=ind) %>%
    select(Sample,total_drops,non_empty,knee_lower=values)


#Write out file. 
write.csv(x = droplet_summary,
          file = here("processed-data","21_Reiner_snRNAseq","snRNA_droplet_summary.csv"),
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

ggsave(plot = droplet_barplot,filename = here("plots","21_Reiner_snRNAseq","droplet_barplot_per_sample.png"))
print("Droplet barplot complete")

#Load in the sce object
sce <- readRDS(here("processed-data","21_Reiner_snRNAseq","sce_raw.Rds"))

dim(sce)

#### Eliminate empty droplets ####
e.out.all <- do.call("rbind", e.out)[colnames(sce), ]
sce <- sce[, which(e.out.all$FDR <= 0.001)]

dim(sce)

#Save object
saveRDS(sce,file = here("processed-data","21_Reiner_snRNAseq","sce_emptyDrops_removed.Rds"))

####Begin QC
sce <- scuttle::addPerCellQC(sce,subsets = list(Mito=which(seqnames(sce) == "MT")))

#Plot mitochondria vs detected
#All samples together
mito_vs_detected <- plotColData(object = sce,
                                y = "subsets_Mito_percent",
                                x = "detected",
                                colour_by = "Sample")

ggsave(filename = here("plots","21_Reiner_snRNAseq","nuclei_QC","mito","mito_vs_detected_bySample.png"),
       plot = mito_vs_detected)
print("mito vs detected plot done")

#Now samples separately.
for(i in unique(sce$Sample)){
    print(i)
    x <- plotColData(object = sce[,sce$Sample == i],
                     y = "subsets_Mito_percent",
                     x = "detected",
                     colour_by = "Sample")
    ggsave(filename = here("plots","21_Reiner_snRNAseq","nuclei_QC","mito",
                           paste0("mito_vs_detected_",i,"_only.png")),
           plot = x)
}

#Before we go on factorize the sample
sce$Sample <- factor(x = sce$Sample)


#########################################
############### High mito ###############
#########################################
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher", batch = sce$Sample)
table(sce$Sample,sce$high_mito)

nMAD3_vln <- plotColData(sce,x = "Sample",y= "subsets_Mito_percent",colour_by = "high_mito") +
    ggtitle(paste0("Mito Percent\n3 MADs")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 5,lty = 2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = nMAD3_vln,filename = here("plots","21_Reiner_snRNAseq",
                                        "nuclei_QC","mito",
                                        "nMAD3_Mito_Pct_violin.png"))

#Simialr to the Comment in https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/51d15ef9f5f2c4c53f55e22e3fe467de1a724668/10x_all-FACS-n10_2021rev_step01_processing-QC_MNT.R#L4
#MAD approach is unnecessarily throwing out cells because the distribution is centered around 0

#Check the number of nuclei 
sce$high_mito_pct <- ifelse(sce$subsets_Mito_percent > 5.0,
                            TRUE,
                            FALSE)

print("Table for number of nuclei removed by high mitochondrial percentage")
table(sce$high_mito_pct)

mito_pct_vln <- plotColData(sce,x = "Sample",y= "subsets_Mito_percent",colour_by = "high_mito_pct") +
    ggtitle(paste0("Mito Percent\n5% cutoff")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 5,lty = 2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = mito_pct_vln,filename = here("plots","21_Reiner_snRNAseq","nuclei_QC","mito","FivePct_Mito_Pct_violin.png"))

#########################################
######### Low detected genes ############
#########################################
# ## low detected genes
#Plot # of detected genes per sample to look at distributions. 
detected_violin_sample <- plotColData(sce, x = "Sample", y = "detected",colour_by = "Sample") +
    scale_y_log10() +
    ggtitle("Number of genes/nucleus") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") 

ggsave(plot = detected_violin_sample,filename = here("plots","21_Reiner_snRNAseq",
                                                     "nuclei_QC","number_genes",
                                                     "number_of_genes_coloredbySample.png"))

sce$low_detected <- isOutlier(sce$detected, nmads = 2, type = "lower", batch = sce$Sample)

low_genes <- plotColData(sce, x = "Sample", y = "detected",colour_by = "low_detected") +
  scale_y_continuous(trans="log10") +
  ggtitle("Number of genes/nucleus\nnMAD=2") +
  geom_hline(yintercept = 750,lty = 2,color = "red") +
  geom_hline(yintercept = 1000,lty = 2,color = "blue") +
  geom_hline(yintercept = 2000,lty = 2,color = "green") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
ggsave(plot = low_genes,here("plots","21_Reiner_snRNAseq",
                               "nuclei_QC","number_genes",
                               "number_genes_nMAD2_batchbySample.png"))

print("Number of nuclei with low genes (nMAD=2)")
table(sce$low_detected,sce$Sample)

print("detected QC plots complete")
#########################################
########## Low Library Size #############
#########################################
# ## low library size
#Plot library size per sample to look at distributions. 
lib_size_violin <- plotColData(sce, x = "Sample", y = "sum",colour_by = "Sample") +
    scale_y_log10() +
    ggtitle("Total UMIs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = lib_size_violin,filename = here("plots","21_Reiner_snRNAseq",
                                              "nuclei_QC","library_size",
                                              "library_size_violin_coloredbySample.png"))


sce$low_library <- isOutlier(sce$sum, nmads = 2, type = "lower", batch = sce$Sample,log = TRUE)

low_lib <- plotColData(sce, x = "Sample", y = "sum",colour_by = "low_library") +
  scale_y_continuous(trans="log2") +
  ggtitle("Total UMIs\nnMAD=2") +
  geom_hline(yintercept = 1000,lty = 2,color = "red") +
  geom_hline(yintercept = 2000,lty = 2,color = "blue") +
  geom_hline(yintercept = 4000,lty = 2,color = "green") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
ggsave(plot = low_lib,here("plots","21_Reiner_snRNAseq",
                             "nuclei_QC","library_size",
                             "library_size_nMAD2_batchbySort_logTrue.png"))

print("Number of nuclei with low library size by sample (nMAD=2)")
table(sce$low_library,sce$Sample)

print("library size QC completed")
########################################################

#### Doublet detection ####
## To speed up, run on sample-level top-HVGs - just take top 1000
print("Running Doublet Detection")
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

ggsave(dbl_box_plot, filename = here("plots","21_Reiner_snRNAseq",
                                     "nuclei_QC","doublet_score",
                                     "doublet_scores_boxplot.png"))

dbl_density_plot <- dbl_df %>%
    ggplot(aes(x = doubletScore,fill = Sample)) +
    geom_density() +
    labs(x = "doublet score") +
    theme_bw()

ggsave(dbl_density_plot, filename = here("plots","21_Reiner_snRNAseq",
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

#Identify cells that were low quality as determined by mitochondrial percentage + above metrics
qc_sum_detected <- sce$sum < 2000 | sce$detected < 750
qc_mito_pct <- sce$subsets_Mito_percent > 5
qc_doublet <- sce$doubletScore >= 5

#Add discard to sce object
sce$discard <- qc_sum_detected | qc_mito_pct | qc_doublet 

print("How many cells are removed from the object?")
table(sce$discard)

qc_t <- addmargins(table(sce$Sample, sce$discard))

qc_t

#Write out qc table
write.csv(qc_t,file = here("processed-data","21_Reiner_snRNAseq","DiscardedCells_by_Sample.csv"))

print("How many cells identified as low quality cells are doublets?")
table(sce$discard,sce$doubletScore >= 5)

#Save object
saveRDS(sce,file = here("processed-data","21_Reiner_snRNAseq","sce_emptyDrops_removed_withQC.Rds"))

#Remove the cells that don't meet the basic QC cutoffs 
sce <- sce[,!sce$discard]
dim(sce)

#Print the sce object. 
print("SCE object with low quality cells and doublets removed")
sce

#Save object
saveRDS(sce,file = here("processed-data","21_Reiner_snRNAseq","sce_clean.Rds"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
