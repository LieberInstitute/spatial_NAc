library(SpatialExperiment)
library(spatialLIBD)
library(jaffelab)
library(here)
library(sessioninfo)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(fastTopics)

dat_dir <- here("processed-data", "12_snRNA")
res_dir <- here("processed-data", "10_post_clustering_analysis", "02_spatial_registration_sn")
plot_dir <- here("plots", "10_post_clustering_analysis", "02_spatial_registration_sn")


#### Load sn data & visualize the number of cells which belong to each type ####
sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
spe_pseudo <- readRDS(file = file.path("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/10_post_clustering_analysis/01_pseudobulk_markers/01_precast/pseudobulk_capture_area/final_clusters/spe_pseudo_precast_clusters.rds"))

geneData <- rowData(spe_pseudo)
geneData <- geneData[match(rowData(sce)$gene_id, geneData$gene_id), ]
rowData(sce)$gene_search <- geneData$gene_search

sce <- sce[ ,!sce$CellType.Final == "Neuron_Ambig"]

subset_neurons <- TRUE
if(subset_neurons){
    sce <- sce[ ,sce$CellType.Final %in% c("DRD1_MSN_A", "DRD1_MSN_B", "DRD1_MSN_C", 
    "DRD1_MSN_D", "DRD2_MSN_A", "DRD2_MSN_B", "Inh_A", "Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F", 
    "Excitatory")]
}

counts <- counts(sce)
cellType_abundance <- data.frame(table(sce$CellType.Final))
colnames(cellType_abundance) <- c("Cell_type", "Ncells")
pdf(file.path(plot_dir, "cellType_abundance.pdf"), width = 12, height = 5)
ggplot(cellType_abundance, aes(x = Cell_type, y = Ncells, fill = Cell_type)) + geom_bar(stat = "identity") + theme_classic() + xlab("Cell Type") + ylab("Number of cells") + coord_flip() + theme(legend.position = "none") 
dev.off()

# Add age and sex information to the sce object
sample_info_path <- here("raw-data", "sample_key_spatial_NAc.csv")
sample_info <- read.csv(sample_info_path) |>
    as_tibble() |>
    select(c(Brain, Age, Sex)) |>
    distinct_all()

Age <- c()
Sex <- c()
for(i in c(1:dim(sce)[2])){
    Age[i] <- sample_info$Age[sample_info$Brain == colData(sce)$Brain_ID[i]]
    Sex[i] <- sample_info$Sex[sample_info$Brain == colData(sce)$Brain_ID[i]]
}

colData(sce)$Sex <- Sex
colData(sce)$Age <- Age

## Factor categorical variables used as covariates
colData(sce)$Sex <- factor(colData(sce)$Sex, levels = c("F", "M"))
sce$age_scaled <- scales::rescale(sce$Age,to=c(0,1))
## Use all unique ensembl IDs as rownames
rownames(sce) <- rowData(sce)$gene_id

# Rename T-cell to avoid downstream issues
sce$CellType.Final[sce$CellType.Final == "T-Cell"] <- "T_cell"

## Run single cell registration to identify markers
sn_registration <- registration_wrapper(
    sce = sce,
    var_registration = "CellType.Final",
    var_sample_id = "Brain_ID",
    covars = c("age_scaled", "Sex"),
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
)

if(subset_neurons){
    saveRDS(sn_registration, file = file.path(res_dir, "sn_cellType_registration_neurons.rds"))
}else{
    saveRDS(sn_registration, file = file.path(res_dir, "sn_cellType_registration.rds"))
}

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()