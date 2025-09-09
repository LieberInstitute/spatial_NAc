####snRNA-seq NMF pattern projection to Visium
library(tidyverse)
library(RcppML)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(here)
library(sessioninfo)
library(scater)
library(scran)
library(BiocParallel)
library(BiocSingular)
library(spatialNAcUtils)
library(jaffelab)
library(projectR)
library(scater)
library(scran)
library(dittoSeq)
library(escheR)
library(getopt)
library(Seurat)

spec <- matrix(c(
  'data',        'd', 1, 'character', 'Specify the input dataset',
  'nFactors',    'n', 1, 'integer',   'Number of NMF factors to use',
  'select_HVGs', 's', 1, 'logical',   'Whether to use highly variable genes (TRUE/FALSE)'
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

#opt <- list()
#opt$data <- "rat_case_control_morphine_acute"
#opt$nFactors <- 30
#opt$select_HVGs <- FALSE

res_dir <- here::here("processed-data", "16_transfer_learning", "04_cellType_NMF")
plot_dir <- here::here("plots", "16_transfer_learning", "04_cellType_NMF")
res_dir <- paste0(res_dir, "/", opt$data)
plot_dir <- paste0(plot_dir, "/", opt$data)

if(opt$select_HVGs){
    x <- readRDS(file.path(res_dir, paste0("nmf_results_HVGs_", opt$nFactors, "factors.rds")))
}else{
    x <- readRDS(file.path(res_dir, paste0("nmf_results_", opt$nFactors, "factors.rds")))
}

patterns <- t(x@h)
colnames(patterns) <- paste("NMF", 1:dim(patterns)[2], sep = "_")
# extract loadings
loadings <- x@w

# Read in and process the human NAc data
raw_in_path <- here(
    "processed-data", "05_harmony_BayesSpace", "02-compute_QC_metrics", "spe_with_QC_metrics_hdf5"
)
spe <- loadHDF5SummarizedExperiment(raw_in_path)
rownames(spe) <- rowData(spe)$gene_name
rownames(spe) <- make.names(rownames(spe), unique = TRUE)
spe$low_umi <- spe$sum_umi < 250
spe$low_gene <- spe$sum_gene < 250
spe$low_gene_edge_spot <- spe$low_umi & spe$edge_distance < 6
spe <- spe[ ,spe$low_gene_edge_spot == "FALSE"]
spe <- spe[ ,spe$local_outliers == "FALSE"]
exprLogic <- counts(spe) > 0
spe$sample_id <- as.factor(spe$sample_id)
select.genes <-  readRDS(file.path(here::here("processed-data", "16_transfer_learning", "04_cellType_NMF"), "select_genes.rds"))
# Only select those genes which have non-zero expression in atleast 1 spot in each slide
spe <- spe[rownames(spe) %in% select.genes, ]

# Now process the orthologs
rat_human_ortholog <- read.delim("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/HOM_AllOrganism.rpt")
rat_orthologs <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "rat", ]
human_orthologs <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "human", ]
common_DB_class_keys <- intersect(rat_orthologs$DB.Class.Key, human_orthologs$DB.Class.Key)
rat_orthologs <- rat_orthologs[rat_orthologs$DB.Class.Key %in% common_DB_class_keys, ]
human_orthologs <- human_orthologs[human_orthologs$DB.Class.Key %in% common_DB_class_keys, ]

ortholog_list <- lapply(common_DB_class_keys, function(id){
    rat_genes <- rat_orthologs$Symbol[rat_orthologs$DB.Class.Key == id]
    rat_genes <- rat_genes[rat_genes %in% rownames(loadings)]
    human_genes <- human_orthologs$Symbol[human_orthologs$DB.Class.Key == id]
    human_genes <- human_genes[human_genes %in% rownames(spe)]
    list("rat_genes" = rat_genes, "human_genes" = human_genes)
})

names(ortholog_list) <- common_DB_class_keys
length_rat_genes <- lapply(ortholog_list, function(iorg){
    length(iorg$rat_genes)
}) 
length_rat_genes <- unlist(length_rat_genes)

length_human_genes <- lapply(ortholog_list, function(iorg){
    length(iorg$human_genes)
}) 
length_human_genes <- unlist(length_human_genes)

ortholog_list <- ortholog_list[!((length_rat_genes == 0) | (length_human_genes == 0))]
length_rat_genes <- lapply(ortholog_list, function(iorg){
    length(iorg$rat_genes)
}) 
length_rat_genes <- unlist(length_rat_genes)

length_human_genes <- lapply(ortholog_list, function(iorg){
    length(iorg$human_genes)
}) 
length_human_genes <- unlist(length_human_genes)
ortholog_list_unique <- ortholog_list[(length_rat_genes == 1) & (length_human_genes == 1)]
rat_genes <- lapply(ortholog_list_unique, function(iorg){
    iorg$rat_genes
})
human_genes <- lapply(ortholog_list_unique, function(iorg){
    iorg$human_genes
})
rat_genes <- unlist(rat_genes)
human_genes <- unlist(human_genes)
orthologs_unique_df <- data.frame("JAX_id" = names(rat_genes), "rat_genes" = rat_genes, "human_genes" = human_genes)
rownames(orthologs_unique_df) <- c(1:dim(orthologs_unique_df)[1])

orthologs_df <- orthologs_unique_df[!duplicated(orthologs_unique_df$rat_genes), ]
orthologs_df <- orthologs_df[!duplicated(orthologs_df$human_genes), ]
orthologs_df <- orthologs_df[orthologs_df$rat_genes %in% rownames(loadings), ]
orthologs_df <- orthologs_df[orthologs_df$human_genes %in% rownames(spe), ]

loadings <- loadings[rownames(loadings) %in% orthologs_df$rat_genes, ]
loadings <- loadings[match(orthologs_df$rat_genes, rownames(loadings)), ]
spe <- spe[rownames(spe) %in% orthologs_df$human_genes, ]
spe <- spe[match(orthologs_df$human_genes, rownames(spe)), ]

# Now get the SPE data
counts <- as.matrix(counts(spe))
metadata <- colData(spe)
sp_obj <- CreateSeuratObject(counts = counts, project = "NAc_10x_visium")
metadata <- metadata[rownames(metadata) %in% colnames(sp_obj), ]
metadata <- metadata[match(colnames(sp_obj), rownames(metadata)), ]
sp_obj$donor <- metadata$sample_id
sp_obj$capture_area <- metadata$sample_id_original
sp_obj$slide_num <- metadata$slide_num
sp_obj$array_num <- metadata$array_num
sp_obj$age <- metadata$Age
sp_obj$sex <- metadata$Sex
sp_obj$exclude_overlapping <- metadata$exclude_overlapping
sp_obj$overlap_slide <- metadata$overlap_slide
sp_obj <- sp_obj %>% 
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000)
logcounts(spe) <- sp_obj[["RNA"]]$data
##projection
set.seed(1029)
# extract patterns
logcounts <- logcounts(spe)

proj <- project(w=loadings, data=logcounts)

# remove rowSums == 0
proj1 <- proj[rowSums(proj) == 0, ,drop = FALSE]
proj2 <- proj[rowSums(proj) != 0, ,drop = FALSE]

proj2 <- apply(proj2,1,function(x){x/sum(x)})
proj1 <- t(proj1)

proj_final <- cbind(proj2, proj1)
proj_final <- proj_final[ ,match(rownames(proj), colnames(proj_final))]

colData(spe) <- cbind(colData(spe),proj_final)
if(opt$select_HVGs){
    saveRDS(spe, file = file.path(res_dir, paste0("spe_NMF_", opt$nFactors ,"_HVGs.rds")))
}else{
    saveRDS(spe, file = file.path(res_dir, paste0("spe_NMF_", opt$nFactors ,".rds")))
}
