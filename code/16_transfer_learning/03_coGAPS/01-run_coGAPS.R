rm(list = ls())
library(Seurat)
library(Matrix)
library(CoGAPS)
library(here)
library(ggplot2)
library(pheatmap)
library(HDF5Array)
library(SpatialExperiment)
library(spatialLIBD)
library(scran)
library(scuttle)
library(SingleCellExperiment)

opt <- list()
opt$data <- "rat_case_control_cocaine"

# Set the locations of the output files and plots
res_dir <- here::here("processed-data", "16_transfer_learning", "03_coGAPS")
plot_dir <- here::here("plots", "16_transfer_learning", "03_coGAPS")
res_dir <- paste0(res_dir, "/", opt$data)
plot_dir <- paste0(plot_dir, "/", opt$data)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

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
nSpots_by_donor <- lapply(levels(spe$sample_id), function(iSample){
    rowSums(exprLogic[ ,spe$sample_id == iSample])
})
nSpots_by_donor <- do.call(cbind, nSpots_by_donor)
colnames(nSpots_by_donor) <- levels(spe$sample_id)
select.genes <-  rownames(nSpots_by_donor)[rowSums(nSpots_by_donor == 0) == 0]
# Only select those genes which have non-zero expression in atleast 1 spot in each slide
spe <- spe[rownames(spe) %in% select.genes, ]
spe <- logNormCounts(spe)

# Read in the single nucleus data that is to be factorized
if(opt$data == "rat_case_control_cocaine_acute"){
    scDir <- here::here("processed-data", "12_snRNA")
    sce <- readRDS(file = file.path(scDir, "NAc_Combo_Integrated.RDS"))
    sce <- subset(sce, Dataset == "Acute")
  }else{
    if(opt$data == "rat_case_control_cocaine_repeated"){
        scDir <- here::here("processed-data", "12_snRNA")
        sce <- readRDS(file = file.path(scDir, "NAc_Combo_Integrated.RDS"))
        sce <- subset(sce, Dataset == "Repeated")
    }else{
      if(opt$data == "rat_case_control_morphine_acute"){
        scDir <- here("processed-data", "21_Reiner_snRNAseq")
        sce <- readRDS(file = file.path(scDir, "sce_clean_w_labels.Rds"))
        sce <- sce[ ,sce$group == "acute"]
      }else{
        if(opt$data == "rat_case_control_morphine_repeated"){
          scDir <- here("processed-data", "21_Reiner_snRNAseq")
          sce <- readRDS(file = file.path(scDir, "sce_clean_w_labels.Rds"))
          sce <- sce[ ,sce$group == "chronic"]
        }else{
          stop("Invalid input data set")
        }
      }
    }  
  }

if(opt$data == "rat_case_control_acute" | opt$data == "rat_case_control_repeated"){
  
  # For dimensionality reduction we only want to use genes that map to human orthologs and are found in the spatial data
  rat_human_ortholog <- read.delim("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/HOM_AllOrganism.rpt")
  rat_orthologs <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "rat", ]
  human_orthologs <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "human", ]
  common_DB_class_keys <- intersect(rat_orthologs$DB.Class.Key, human_orthologs$DB.Class.Key)
  rat_orthologs <- rat_orthologs[rat_orthologs$DB.Class.Key %in% common_DB_class_keys, ]
  human_orthologs <- human_orthologs[human_orthologs$DB.Class.Key %in% common_DB_class_keys, ]
  ortholog_list <- lapply(common_DB_class_keys, function(id){
    rat_genes <- rat_orthologs$Symbol[rat_orthologs$DB.Class.Key == id]
    rat_genes <- rat_genes[rat_genes %in% rownames(sce)]
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
  orthologs_df <- orthologs_df[orthologs_df$rat_genes %in% rownames(sce), ]
  orthologs_df <- orthologs_df[orthologs_df$human_genes %in% rownames(spe), ]
  saveRDS(orthologs_df, file.path(res_dir, "orthologs_df.rds"))

  # Read the differentially expressed genes in cocaine vs. saline
  diffExpr_file <- here("processed-data", "16_transfer_learning", "rat_case_control_DEGs", "NAc_rat_cocaine_all.csv")
  DEGs_cocaine <- read.csv(diffExpr_file)
  print(length(intersect(unique(DEGs_cocaine$Gene), orthologs_df$rat_genes)))

  sce <- sce[rownames(sce) %in% orthologs_df$rat_genes, ]
  
  subsample <- FALSE
  if(subsample){
    meta <- sce@meta.data
    meta$group <- paste(meta$Stim, meta$Combo_CellType, sep = "_")
    group_sizes <- table(meta$group)
    min_cells <- 200
    set.seed(42)  # for reproducibility
    sampled_cells <- unlist(lapply(names(group_sizes), function(g) {
    group_cells <- rownames(meta)[meta$group == g]
    if (length(group_cells) >= min_cells) {
        sample(group_cells, min_cells)
    } else {
        group_cells  # if group smaller than min_cells, include all (optional)
    }
    }))
    sce <- subset(sce, cells = sampled_cells)
  }  

  sce <- NormalizeData(sce)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, features = VariableFeatures(sce))
  
  var_genes <- VariableFeatures(sce)
  print(length(intersect(unique(DEGs_cocaine$Gene), var_genes)))

  data_matrix <- GetAssayData(sce, assay = "RNA", slot = "data")[var_genes, ]

  # --- PCA ---
  sce <- RunPCA(sce, features = VariableFeatures(sce))

  # --- UMAP ---
  sce <- RunUMAP(sce, dims = 1:20)

  pdf(file.path(plot_dir, "UMAP_cellType.pdf"), width = 6, height = 4)
  DimPlot(sce, group.by = "Combo_CellType")
  dev.off()

  pdf(file.path(plot_dir, "UMAP_Stim.pdf"), width = 6, height = 4)
  DimPlot(sce, group.by = "Stim")
  dev.off()

  pdf(file.path(plot_dir, "UMAP_Dataset.pdf"), width = 6, height = 4)
  DimPlot(sce, group.by = "Dataset")
  dev.off()
  
  saveRDS(sce, file.path(res_dir, "snRNA_seq.rds"))
}

if(opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
  rownames(sce) <- rowData(sce)$gene_name
  rownames(sce) <- make.names(rownames(sce), unique = TRUE)
  counts <- counts(sce)
  metadata <- colData(sce)
  sobj <- CreateSeuratObject(counts = counts, project = "rat_morphine", min.cells = 1)
  metadata <- metadata[rownames(metadata) %in% colnames(sobj), ]
  metadata <- metadata[match(colnames(sobj), rownames(metadata)), ]
  sobj$Sample <- metadata$Sample
  sobj$sum <- metadata$sum
  sobj$detected <- metadata$detected
  sobj$percent_mito <- metadata$subsets_Mito_percent
  sobj$CellType <- metadata$cluster
  sobj$group <- metadata$group
  sobj$treatment <- metadata$treatment
  sobj$subject <- metadata$subject
  sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)
 
  # For dimensionality reduction we only want to use genes that map to human orthologs and are found in the spatial data
  rat_human_ortholog <- read.delim("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/HOM_AllOrganism.rpt")
  rat_orthologs <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "rat", ]
  human_orthologs <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "human", ]
  common_DB_class_keys <- intersect(rat_orthologs$DB.Class.Key, human_orthologs$DB.Class.Key)
  rat_orthologs <- rat_orthologs[rat_orthologs$DB.Class.Key %in% common_DB_class_keys, ]
  human_orthologs <- human_orthologs[human_orthologs$DB.Class.Key %in% common_DB_class_keys, ]
  ortholog_list <- lapply(common_DB_class_keys, function(id){
    rat_genes <- rat_orthologs$Symbol[rat_orthologs$DB.Class.Key == id]
    rat_genes <- rat_genes[rat_genes %in% rownames(sobj)]
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
  orthologs_df <- orthologs_df[orthologs_df$rat_genes %in% rownames(sobj), ]
  orthologs_df <- orthologs_df[orthologs_df$human_genes %in% rownames(spe), ]
  saveRDS(orthologs_df, file.path(res_dir, "orthologs_df.rds"))

  # Read the differentially expressed genes in cocaine vs. saline
  if(opt$data == "rat_case_control_morphine_acute"){
    diffExpr_file <- here("processed-data", "21_Reiner_snRNAseq", "DEGs_acute.csv")
  }
  if(opt$data == "rat_case_control_morphine_repeated"){
    diffExpr_file <- here("processed-data", "21_Reiner_snRNAseq", "DEGs_chronic.csv")
  }
  DEGs_morphine <- read.csv(diffExpr_file)
  print("Number of unique DEGs:")
  print(length(unique(DEGs_morphine$Gene)))
  print("Number of DEGs included after selecting orthologs:")
  print(length(intersect(unique(DEGs_morphine$Gene), orthologs_df$rat_genes)))


  sobj <- sobj[rownames(sobj) %in% orthologs_df$rat_genes, ]
  # Select genes for which we have orthologs
  sobj <- sobj %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)
}


##############################################################################
# Run coGAPS
##############################################################################
params <- CogapsParams(nIterations=10000, # run for 10000 iterations 
               	       seed=42, # for consistency across stochastic runs
               	       nPatterns=30, # each thread will learn 8 patterns
                       sparseOptimization=TRUE, # optimize for sparse data
                       distributed="genome-wide") # parallelize across sets
params <- setDistributedParams(params, nSets=10)
cogaps_result <- CoGAPS(data = as.matrix(data_matrix), params = params)

saveRDS(cogaps_result, file.path(res_dir, "cogaps_results.rds"))