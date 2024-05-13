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
library(readr)
library(scran)
library(scater)
library(purrr)
library(dplyr)
library(tidyr)
library(here)

#Compile droplet qc information
print("Reading in droplet paths")
droplet_paths <- list.files(path = here("processed-data","12_snRNA","droplet_scores"),
                            full.names = TRUE)

names(droplet_paths) <- gsub(x = basename(droplet_paths),
                             pattern = "_droplet_scores.Rdata",
                             replacement = "")

droplet_paths


#Read in the droplet scores
e.out <- lapply(droplet_paths, function(x) get(load(x)))

#To make sure we aren't throwing out any cells check if Limited=TRUE and SIG==FALSE
#If both are true, then we could be throwing out non-empty droplets. 
map(e.out, ~ addmargins(table(Signif = .x$FDR <= 0.001, Limited = .x$Limited)))


#Pull knee lower values
std_out <- list.files(here("code", "12_build_sce", "logs"), pattern = "emptydrops", full.names = TRUE)
std_out <- map(std_out, readLines)

#160th line of each iteration of the list is the knee_lower value
knee_lower <- as.character(lapply(X = std_out,"[",160))
knee_lower

#Just keep the values
knee_lower <- as.numeric(lapply(strsplit(knee_lower,split = "="),"[",2))


#151st element of the std_out list contains the sample na,e
sample_name <- as.character(lapply(std_out,"[",151))
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


droplet_summary

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
print("Droplet barplot complete")

#Load in the sce object
sce <- readRDS(here("processed-data","12_snRNA","sce_raw.Rds"))

dim(sce)

#### Eliminate empty droplets ####
e.out.all <- do.call("rbind", e.out)[colnames(sce), ]
sce <- sce[, which(e.out.all$FDR <= 0.001)]

dim(sce)

#Save object
saveRDS(sce,file = here("processed-data","12_snRNA","sce_emptyDrops_removed.Rds"))

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
print("mito vs detected plot done")

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

#Before we go on factorize the sample
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
ggsave(plot = nMAD3_vln,filename = here("plots","12_snRNA",
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
ggsave(plot = mito_pct_vln,filename = here("plots","12_snRNA","nuclei_QC","mito","FivePct_Mito_Pct_violin.png"))

#########################################
######### Low detected genes ############
#########################################
# ## low detected genes
#Plot # of detected genes per sample to look at distributions. 
detected_violin_sample <- plotColData(sce, x = "Sample", y = "detected",colour_by = "Sample") +
    scale_y_log10() +
    ggtitle("Number of genes/nucleus") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = detected_violin_sample,filename = here("plots","12_snRNA",
                                                     "nuclei_QC","number_genes",
                                                     "number_of_genes_coloredbySample.png"))

#Plot # of detected genes and color by sort
detected_violin_sort <- plotColData(sce, x = "Sample", y = "detected",colour_by = "Sort") +
    scale_y_log10() +
    ggtitle("Number of genes/nucleus") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = detected_violin_sort,filename = here("plots","12_snRNA",
                                                   "nuclei_QC","number_genes",
                                                   "number_of_genes_coloredbySort.png"))

#Plot # of detected genes and color by Brain
detected_violin_brainID <- plotColData(sce, x = "Sample", y = "detected",colour_by = "Brain_ID") +
    scale_y_log10() +
    ggtitle("Number of genes/nucleus") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = detected_violin_brainID,filename = here("plots","12_snRNA",
                                                      "nuclei_QC","number_genes",
                                                      "number_of_genes_coloredbyBrainID.png"))

#NeuN sorted samples are unimodal, while PI sorted samples are bimodal. 
#Will batch by Sort type because of this. 
#MAD approach with batch by Sample results in the removal of high quality cells from 
#NeuN sorted samples and retention of low quality cells from PI sorted samples. 
#nMAD=3
sce$low_detected_3 <- isOutlier(sce$detected, nmads = 3, type = "lower", batch = sce$Sort)

low_genes_3 <- plotColData(sce, x = "Sample", y = "detected",colour_by = "low_detected_3") +
  scale_y_continuous(trans="log10") +
  ggtitle("Number of genes/nucleus\nnMAD=3") +
  geom_hline(yintercept = 750,lty = 2,color = "red") +
  geom_hline(yintercept = 1000,lty = 2,color = "blue") +
  geom_hline(yintercept = 2000,lty = 2,color = "green") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
ggsave(plot = low_genes_3,here("plots","12_snRNA",
                               "nuclei_QC","number_genes",
                               "number_genes_nMAD3_batchbySort.png"))

print("Number of nuclei with low genes by sort (nMAD=3)")
table(sce$low_detected_3,sce$Sort)

#nMAD=2
sce$low_detected_2 <- isOutlier(sce$detected, nmads = 2, type = "lower", batch = sce$Sort)

low_genes_2 <- plotColData(sce, x = "Sample", y = "detected",colour_by = "low_detected_2") +
  scale_y_continuous(trans="log10") +
  ggtitle("Number of genes/nucleus\nnMAD=2") +
  geom_hline(yintercept = 750,lty = 2,color = "red") +
  geom_hline(yintercept = 1000,lty = 2,color = "blue") +
  geom_hline(yintercept = 2000,lty = 2,color = "green") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
ggsave(plot = low_genes_2,here("plots","12_snRNA",
                               "nuclei_QC","number_genes",
                               "number_genes_nMAD2_batchbySort.png"))

print("Number of nuclei with low genes by sort (nMAD=2)")
table(sce$low_detected_2,sce$Sort)

#nMAD=1
sce$low_detected_1 <- isOutlier(sce$detected, nmads = 1, type = "lower", batch = sce$Sort)

low_genes_1 <-plotColData(sce, x = "Sample", y = "detected",colour_by = "low_detected_1") +
  scale_y_continuous(trans="log10") +
  ggtitle("Number of genes/nucleus\nnMAD=1") +
  geom_hline(yintercept = 1000,lty = 2,color = "red") +
  geom_hline(yintercept = 2000,lty = 2,color = "blue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
ggsave(plot = low_genes_1,here("plots","12_snRNA",
                               "nuclei_QC","number_genes",
                               "number_genes_nMAD1_batchbySort.png"))

print("Number of nuclei with low genes by sort (nMAD=1)")
table(sce$low_detected_1,sce$Sort)

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
ggsave(plot = lib_size_violin,filename = here("plots","12_snRNA",
                                              "nuclei_QC","library_size",
                                              "library_size_violin_coloredbySample.png"))

#Plot library size per sort 
lib_size_violin_sort <- plotColData(sce, x = "Sample", y = "sum",colour_by = "Sort") +
    scale_y_log10() +
    ggtitle("Total UMIs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = lib_size_violin_sort,filename = here("plots","12_snRNA",
                                              "nuclei_QC","library_size",
                                              "library_size_violin_coloredbySort.png"))

#Plot library size per brain 
lib_size_violin_brain <- plotColData(sce, x = "Sample", y = "sum",colour_by = "Brain_ID") +
    scale_y_log10() +
    ggtitle("Total UMIs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = lib_size_violin_brain,filename = here("plots","12_snRNA",
                                                    "nuclei_QC","library_size",
                                                    "library_size_violin_coloredbyBrain.png"))

#Batching by Sort again for library size. Also, log=TRUE so library size is log2 transformed before
#MADs are computed.
#nMAD=3
sce$low_library_3 <- isOutlier(sce$sum, nmads = 3, type = "lower", batch = sce$Sort,log = TRUE)

low_lib_3 <- plotColData(sce, x = "Sample", y = "sum",colour_by = "low_library_3") +
  scale_y_continuous(trans="log2") +
  ggtitle("Total UMIs\nnMAD=3") +
  geom_hline(yintercept = 1000,lty = 2,color = "red") +
  geom_hline(yintercept = 2000,lty = 2,color = "blue") +
  geom_hline(yintercept = 4000,lty = 2,color = "green") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
ggsave(plot = low_lib_3,here("plots","12_snRNA",
                             "nuclei_QC","library_size",
                             "library_size_nMAD3_batchbySort_logTrue.png"))

print("Number of nuclei with low library size by sort (nMAD=3)")
table(sce$low_library_3,sce$Sort)

#nMAD=2
sce$low_library_2 <- isOutlier(sce$sum, nmads = 2, type = "lower", batch = sce$Sort,log = TRUE)

low_lib_2 <-plotColData(sce, x = "Sample", y = "sum",colour_by = "low_library_2") +
  scale_y_continuous(trans="log2") +
  ggtitle("Total UMIs\nnMAD=2") +
  geom_hline(yintercept = 1000,lty = 2,color = "red") +
  geom_hline(yintercept = 2000,lty = 2,color = "blue") +
  geom_hline(yintercept = 4000,lty = 2,color = "green") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
ggsave(plot = low_lib_2,here("plots","12_snRNA",
                             "nuclei_QC","library_size",
                             "library_size_nMAD2_batchbySort_logTrue.png"))

print("Number of nuclei with low library size by sort (nMAD=2)")
table(sce$low_library_2,sce$Sort)

#nMAD=1
sce$low_library_1 <- isOutlier(sce$sum, nmads = 1, type = "lower", batch = sce$Sort,log = TRUE)

low_lib_1 <- plotColData(sce, x = "Sample", y = "sum",colour_by = "low_library_1") +
  scale_y_continuous(trans="log2") +
  ggtitle("Total UMIs\nnMAD=1") +
  geom_hline(yintercept = 1000,lty = 2,color = "red") +
  geom_hline(yintercept = 2000,lty = 2,color = "blue") +
  geom_hline(yintercept = 4000,lty = 2,color = "green") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
ggsave(plot = low_lib_1,here("plots","12_snRNA",
                             "nuclei_QC","library_size",
                             "library_size_nMAD1_batchbySort_logTrue.png"))

print("Number of nuclei with low library size by sort (nMAD=1)")
table(sce$low_detected_1,sce$Sort) 

print("library size QC completed")
########################################################
###For NeuN+, low detected genes nMAD = 3, low library = 3
###For PI+, low detected genes nMAD = 1, low library = 1
####Identify NeuN cells that need to be dropped 
NeuN_Drop <- sce[,sce$Sort == "PI_NeuN"]$low_detected_3 | sce[,sce$Sort == "PI_NeuN"]$low_library_3
print("How many NeuN sorted cells are dropped?")
table(NeuN_Drop)

#Get cell IDs that need to be dropped
NeuN_Drop_names <- rownames(colData(sce[,sce$Sort == "PI_NeuN"])[which(NeuN_Drop),])

print("How many NeuN sorted cells removed with low detected genes?")
table(colData(sce[,sce$Sort == "PI_NeuN"])[NeuN_Drop_names,"low_detected_3"])

print("How many NeuN sorted cells removed with low library size?")
table(colData(sce[,sce$Sort == "PI_NeuN"])[NeuN_Drop_names,"low_library_3"])

print("How many PI_NeuN nuclei that have low detected genes also have low library size?")
table(colData(sce[,sce$Sort == "PI_NeuN"])[NeuN_Drop_names,"low_library_3"],
      colData(sce[,sce$Sort == "PI_NeuN"])[NeuN_Drop_names,"low_detected_3"])

####Identify PI cells that need to be dropped 
PI_Drop <- sce[,sce$Sort == "PI"]$low_detected_1 | sce[,sce$Sort == "PI"]$low_library_1
print("How many PI sorted cells are dropped?")
table(PI_Drop)

#Get cell IDs that need to be dropped
PI_Drop_names <- rownames(colData(sce[,sce$Sort == "PI"])[which(PI_Drop),])

print("How many PI sorted cells removed with low detected features?")
table(colData(sce[,sce$Sort == "PI"])[PI_Drop_names,"low_detected_1"])

print("How many PI sorted cells removed with low library size?")
table(colData(sce[,sce$Sort == "PI"])[PI_Drop_names,"low_library_1"])

print("How many	PI_NeuN nuclei that have low detected genes also have low library size?")
table(colData(sce[,sce$Sort == "PI"])[PI_Drop_names,"low_library_1"],
      colData(sce[,sce$Sort == "PI"])[PI_Drop_names,"low_detected_1"])

#Add a unique rowname to the sce object to identify cells to drop
sce$unique_rowname <- rownames(colData(sce))
stopifnot(identical(rownames(colData(sce)),sce$unique_rowname))

#add colDAta column to remove the low quality cells as determined by library size and detected genes
cells_to_drop <- c(NeuN_Drop_names,PI_Drop_names)
sce$discard_sum_detected <- ifelse(sce$unique_rowname %in% cells_to_drop,
                                      TRUE,
                                      FALSE)

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

#Identify cells that were low quality as determined by mitochondrial percentage + above metrics
qc_sum_detected <- sce$discard_sum_detected
qc_mito_pct <- sce$subsets_Mito_percent > 5
qc_doublet <- sce$doubletScore >= 5

#Previous clustering analysis identified low-quality neurons that express far fewer genes and exhbiti
#significant enrichment of immediate early genes such as Fos and NPAS4. I saved the IDs of these cells 
#and will remove them as part of the QC. 
load(here("processed-data","12_snRNA","IEG_neurons_to_remove.rda"),verbose = TRUE)

head(IEG_cells)

#add unique rowname column to the colData
colData(sce)$unique_rowname <- rownames(colData(sce))
#Add TRUE/FALSE statement for whether cell is in the IEG_cells vector
sce$IEG_cells <- colData(sce)$unique_rowname %in% IEG_cells
#sanity check for above code
all(subset(colData(sce),subset=(IEG_cells == TRUE))$unique_rowname %in% IEG_cells)
qc_IEG <- sce$IEG_cells

print("Number of cells identified as low quality neurons (IEG cells)")
table(sce$IEG_cells)

#How many cells removed by sum/detected cutoffs?
print("Number of cells removed by sum/detected cutoffs")
table(sce$discard_sum_detected)

#Add discard to sce object
sce$discard <- qc_sum_detected | qc_mito_pct | qc_doublet | qc_IEG

print("How many cells are removed from the object?")
table(sce$discard)

qc_t <- addmargins(table(sce$Sample, sce$discard))

qc_t

#Write out qc table
write.csv(qc_t,file = here("processed-data","12_snRNA","DiscardedCells_by_Sample.csv"))

print("How many cells identified as low quality cells are doublets?")
table(sce$discard,sce$doubletScore >= 5)

#Save object
saveRDS(sce,file = here("processed-data","12_snRNA","sce_emptyDrops_removed_withQC.Rds"))

#Remove the cells that don't meet the basic QC cutoffs 
sce <- sce[,!sce$discard]
dim(sce)

#Print the sce object. 
print("SCE object with low quality cells and doublets removed")
sce

#Save object
saveRDS(sce,file = here("processed-data","12_snRNA","sce_clean.Rds"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
