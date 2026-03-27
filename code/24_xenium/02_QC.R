#explore quality control metrics
# cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
# module load r_nac

library(SpatialExperiment)
library(SpotSweeper) #For plotting functions. 
library(sessioninfo)
library(ggplot2)
library(scuttle)
library(escheR)
library(scater)
library(dplyr)
library(scran)
library(here)

#Read in the raw spe file
spe <- readRDS(here("processed-data","xenium","SPEs","spe_raw.Rds"))
spe

#Find genes for scuttle subsets 
is_neg <- stringr::str_detect(rownames(spe), "^NegControlProbe")
is_neg2 <- stringr::str_detect(rownames(spe), "^NegControlCodeword")
is_unassigned <- stringr::str_detect(rownames(spe), "^Unassigned")
is_anyneg <- is_neg | is_neg2 | is_unassigned
is_GEX <- rowData(spe)$Type == "Gene Expression" #This will help identify QC based on just the genes in panel

#identify a set of genes to classify cells as neuronal or non_neuronal
neuron_genes <- c("GAD1","GAD2","SLC32A1","RBFOX3","SNAP25",
                  "PPP1R1B","SLC17A6","SLC17A7")
non_neuron_genes <- c("AQP4","GJA1","SLC1A2","MOBP",
                      "ST18","CSPG4","PDGFRA","FOXJ1",
                      "ARHGAP15")

is_neuron <- rownames(spe) %in% neuron_genes
is_nonneuron <- rownames(spe) %in% non_neuron_genes

#Per: https://genomics.uci.edu/wp-content/uploads/sites/30/PDF_UCI_GRT_Hub_Xenium_Workshop_20240118-120fec8a1ce5134f.pdf
# "Negative Control Codewords are
# a random subset of codewords,
# with identical properties to gene
# codewords" 
#and
# "By definition, a call made to negative control codeword is an error"
#Unassigned probes help measure noise and off-target activity because they do not match any gene codewords
#Negative control probes are sequences that should not bind anything and therefore measure false positive rates and specificity

#Add QC metrics 
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(negProbe = is_neg,
                                                        negCodeword = is_neg2,
                                                        unassigned = is_unassigned,
                                                        any_neg = is_anyneg,
                                                        GEX = is_GEX,
                                                        neuronal = is_neuron,
                                                        non_neuron = is_nonneuron))


dim(spe)

#Remove cells that are 0 for everything. 
spe <- spe[, colSums(counts(spe)) > 0]

dim(spe)

#Remove any cell with any neg percent > 25
spe <- spe[,spe$subsets_any_neg_percent <= 25]

dim(spe)

#Now remove cells that are 0 for all gene expression
spe <- spe[,spe$subsets_GEX_percent > 0]

dim(spe)

#No cells removed. 
spe$Neuron <- ifelse(spe$subsets_neuronal_sum >0,
                     "Neuron",
                     "Non_Neuron")

table(spe$Neuron)


neuron_cols <- c("dodgerblue","purple")
names(neuron_cols) <- c("Neuron","Non_Neuron")

for(i in unique(spe$Sample)){
  print(i)
  spe_sub <- spe[,spe$Sample == i]
  #Plot total_counts on top of the tissue to visualize rotation. 
  x <- escheR::make_escheR(spe_sub) |>
    escheR::add_fill("Neuron") +
    scale_fill_manual(values = neuron_cols)
  ggsave(plot = x,
         filename = here("plots","xenium","01_buildSPE","QC",
                         paste0(i,"Neuron_NonNeuron.png")),
         height = 16, width = 18)
  
}



#Calculate log10 for detected genes and sum
spe$log10_detected <- log10(spe$detected)
spe$log10_sum <- log10(spe$sum)

#SPlit object by neuron and non-neuron
spe_neuron <- spe[,spe$Neuron == "Neuron"]
spe_nonneuron <- spe[,spe$Neuron != "Neuron"]


# To extract the actual MAD cutoff 
get_values <- function(x, type = "lower",nmad) {
  center <- median(x, na.rm = TRUE)
  spread <- mad(x, na.rm = TRUE)
  cutoff <- center - nmad * spread
  cutoff
}


#Calcualte cutoffs 
#Detected
neuron_cutoffs_detected_4 <- get_values(spe_neuron$log10_detected, nmad = 4)
neuron_cutoffs_detected_3 <- get_values(spe_neuron$log10_detected, nmad = 3)

non_neuron_cutoffs_detected_4 <- get_values(spe_nonneuron$log10_detected,nmad = 4)
non_neuron_cutoffs_detected_3 <- get_values(spe_nonneuron$log10_detected,nmad = 3)

#Sum
neuron_cutoffs_sum_4 <- get_values(spe_neuron$log10_sum, nmad = 4)
neuron_cutoffs_sum_3 <- get_values(spe_neuron$log10_sum, nmad = 3)

non_neuron_cutoffs_sum_4 <- get_values(spe_nonneuron$log10_sum, nmad = 4)
non_neuron_cutoffs_sum_3 <- get_values(spe_nonneuron$log10_sum, nmad = 3)

#CELL AREA
# To extract the actual MAD cutoff 
get_values_cellarea <- function(x, type = "higher",nmad) {
  center <- median(x, na.rm = TRUE)
  spread <- mad(x, na.rm = TRUE)
  cutoff <- center + nmad * spread
  cutoff
}

neuron_cutoffs_cellarea_4 <- get_values_cellarea(spe_neuron$cell_area, nmad = 4)
neuron_cutoffs_cellarea_3 <- get_values_cellarea(spe_neuron$cell_area, nmad = 3)

non_neuron_cutoffs_cellarea_4 <- get_values_cellarea(spe_nonneuron$cell_area, nmad = 4)
non_neuron_cutoffs_cellarea_3 <- get_values_cellarea(spe_nonneuron$cell_area, nmad = 3)



#Neuron histograms
neuron_detected_histogram <- ggplot(colData(spe_neuron), 
                                    aes(x = log10_detected)) + 
  geom_histogram(bins = 50,colour = "dodgerblue") + 
  geom_vline(xintercept = neuron_cutoffs_detected_4,lty = 2, color = "black") +
  geom_vline(xintercept = neuron_cutoffs_detected_3,lty = 2, color = "red") +
  ggtitle(paste0("Neuronal\nGenes Detected\n",
                 "4 MADs = ",round(10^neuron_cutoffs_detected_4,digits = 4),"\n",
                 "3 MADs = ",round(10^neuron_cutoffs_detected_3,digits = 4))) +
  labs(x = "log10(detected)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


neuron_sum_histogram <- ggplot(colData(spe_neuron), 
                               aes(x = log10_sum)) + 
  geom_histogram(bins = 50,colour = "dodgerblue") + 
  geom_vline(xintercept = neuron_cutoffs_sum_4,lty = 2, color = "black") +
  geom_vline(xintercept = neuron_cutoffs_sum_3,lty = 2, color = "red") +
  ggtitle(paste0("Neuronal\nLibrary Size\n",
                 "4 MADs = ",round(10^neuron_cutoffs_sum_4,digits = 4),"\n",
                 "3 MADs = ",round(10^neuron_cutoffs_sum_3,digits = 4))) +
  labs(x = "log10(detected)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


neuron_cell_area_histogram <-  ggplot(colData(spe_neuron), 
                                      aes(x = cell_area)) + 
  geom_histogram(bins = 50,colour = "dodgerblue") + 
  geom_vline(xintercept = neuron_cutoffs_cellarea_4,lty = 2, color = "black") +
  geom_vline(xintercept = neuron_cutoffs_cellarea_3,lty = 2, color = "red") +
  ggtitle(paste0("Neuronal\nCell area\n",
                 "4 MADs = ",round(neuron_cutoffs_cellarea_4,digits = 4),"\n",
                 "3 MADs = ",round(neuron_cutoffs_cellarea_3,digits = 4))) +
  labs(x = "cell_area") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))



#Non-neuron histograms
non_neuron_detected_histogram <- ggplot(colData(spe_nonneuron), 
                                        aes(x = log10_detected)) + 
  geom_histogram(bins = 50,colour = "purple") + 
  geom_vline(xintercept = non_neuron_cutoffs_detected_4,lty = 2, color = "black") +
  geom_vline(xintercept = non_neuron_cutoffs_detected_3,lty = 2, color = "red") +
  ggtitle(paste0("Non-neuronal\nLibrary Size\n",
                 "4 MADs = ",round(10^non_neuron_cutoffs_detected_4,digits = 4),"\n",
                 "3 MADs = ",round(10^non_neuron_cutoffs_detected_3,digits = 4))) +
  labs(x = "log10(detected)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


non_neuron_sum_histogram <- ggplot(colData(spe_nonneuron), 
                                   aes(x = log10_sum)) + 
  geom_histogram(bins = 50,colour = "purple") + 
  geom_vline(xintercept = non_neuron_cutoffs_sum_4,lty = 2, color = "black") +
  geom_vline(xintercept = non_neuron_cutoffs_sum_3,lty = 2, color = "red") +
  ggtitle(paste0("Non-neuronal\nLibrary Size\n",
                 "4 MADs = ",round(10^non_neuron_cutoffs_sum_4,digits = 4),"\n",
                 "3 MADs = ",round(10^non_neuron_cutoffs_sum_3,digits = 4))) +
  labs(x = "log10(sum)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

nonneuron_cell_area_histogram <-  ggplot(colData(spe_nonneuron), 
                                         aes(x = cell_area)) + 
  geom_histogram(bins = 50,colour = "purple") + 
  geom_vline(xintercept = neuron_cutoffs_cellarea_4,lty = 2, color = "black") +
  geom_vline(xintercept = neuron_cutoffs_cellarea_3,lty = 2, color = "red") +
  ggtitle(paste0("Neuronal\nCell area\n",
                 "4 MADs = ",round(non_neuron_cutoffs_cellarea_4,digits = 4),"\n",
                 "3 MADs = ",round(non_neuron_cutoffs_cellarea_3,digits = 4))) +
  labs(x = "cell_area") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))



histos <- cowplot::plot_grid(plotlist = list(neuron_detected_histogram,neuron_sum_histogram,neuron_cell_area_histogram,
                                             non_neuron_detected_histogram,non_neuron_sum_histogram,nonneuron_cell_area_histogram),
                             ncol = 3)
ggsave(histos,filename = here("plots","xenium","01_buildSPE","QC","split_Global_QC_histogram.pdf"),
       height = 12, width = 12)


######Neuron
#nMAD = 4
#detected
spe_neuron$detected_outlier <- isOutlier(
  spe_neuron$detected,
  log = TRUE,
  nmads = 4,
  type = "lower"
)

#sum
spe_neuron$sum_outlier <- isOutlier(
  spe_neuron$sum,
  log = TRUE,
  nmads = 4,
  type = "lower"
)

#cell area
spe_neuron$cell_area_outlier <- isOutlier(
  spe_neuron$cell_area,
  log = FALSE,
  nmads = 3,
  type = "higher"
)


#nmad=3
spe_nonneuron$detected_outlier <- isOutlier(
  spe_nonneuron$detected,
  log = TRUE,
  nmads = 3,
  type = "lower"
)

spe_nonneuron$sum_outlier <- isOutlier(
  spe_nonneuron$sum,
  log = TRUE,
  nmads = 3,
  type = "lower"
)


#cell area
spe_nonneuron$cell_area_outlier <- isOutlier(
  spe_nonneuron$cell_area,
  log = FALSE,
  nmads = 3,
  type = "higher"
)


#Add outlier information to the object before saving
#First pull cells that are neurons or not
neurons <- intersect(colnames(spe),colnames(spe_neuron))
nonneurons <- intersect(colnames(spe),colnames(spe_nonneuron))
stopifnot(sum(length(neurons),length(nonneurons)) == length(colnames(spe)))

#Add detected outlier information
colData(spe)[neurons,"detected_outlier"] <- colData(spe_neuron)[neurons,"detected_outlier"]
colData(spe)[nonneurons,"detected_outlier"] <- colData(spe_nonneuron)[nonneurons,"detected_outlier"]

#Add sum outlier information
colData(spe)[neurons,"sum_outlier"] <- colData(spe_neuron)[neurons,"sum_outlier"]
colData(spe)[nonneurons,"sum_outlier"] <- colData(spe_nonneuron)[nonneurons,"sum_outlier"]

#Add cell area outlier information
colData(spe)[neurons,"cell_area_outlier"] <- colData(spe_neuron)[neurons,"cell_area_outlier"]
colData(spe)[nonneurons,"cell_area_outlier"] <- colData(spe_nonneuron)[nonneurons,"cell_area_outlier"]

#Save spe with QC information 
message(paste0("Saving spe object with QC - ",Sys.time()))
saveRDS(spe,here("processed-data","xenium","SPEs","spe_splitQC_Global_prediscard.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
