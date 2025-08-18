library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(scater)
library(scran)
library(here)
library(scuttle)
library(Matrix)
library(SingleCellExperiment)
library(HDF5Array)
library(SpatialExperiment)
library(spatialLIBD)
library(RcppML)
library(singlet)
library(getopt)
set.seed(1234)

spec <- matrix(
    c("data", "d", 1, "character", "Specify the input dataset"),
    byrow = TRUE, ncol = 5
)

opt <- getopt(spec)

#opt <- list()
#opt$data <- "rat_case_control_repeated"

# Read data and create Seurat object
dat_dir <- here::here("processed-data", "12_snRNA")
res_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "preliminary_analysis")
plot_dir <- here::here("plots", "16_transfer_learning", "01_process_reference", "preliminary_analysis")
res_dir <- paste0(res_dir, "/", opt$data)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
# Read in the spatial data
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

if(opt$data == "human_NAc"){
  sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
}else{
  if(opt$data == "rat_case_control_acute"){
    sce <- readRDS(file = file.path(dat_dir, "NAc_Combo_Acute.RDS"))
  }else{
    if(opt$data == "rat_case_control_repeated"){
      sce <- readRDS(file = file.path(dat_dir, "NAc_Combo_Repeated.RDS"))
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
}

if(opt$data == "human_NAc"){
  # Initialize the Seurat object
  sce <- sce[ ,!sce$CellType.Final == "Neuron_Ambig"]
  counts <- counts(sce)
  metadata <- colData(sce)
  sobj <- CreateSeuratObject(counts = counts, project = "NAc_snRNAseq")
  metadata <- metadata[rownames(metadata) %in% colnames(sobj), ]
  metadata <- metadata[match(colnames(sobj), rownames(metadata)), ]
  sobj$Sample <- metadata$Sample
  sobj$Brain_ID <- metadata$Brain_ID
  sobj$percent_mito <- metadata$subsets_Mito_percent
  sobj$CellType <- metadata$CellType.Final
  sobj$Sort <- metadata$Sort
  sobj$snRNA_date <- metadata$snRNA_data
  sobj$Chromium_cDNA_date <- metadata$Chromium_cDNA_date
  sobj$Chromium_library_date <- metadata$Chromium_library_date
  sobj$Chromium_kit <- metadata$Chromium_kit
  sobj$Library_Kit <- metadata$Library_Kit

  # Only select genes that were expressed in the spatial data
  exprLogic2 <- sobj[["RNA"]]$counts > 0
  sobj$Brain_ID <- as.factor(sobj$Brain_ID)
  nCells_by_donor <- lapply(levels(sobj$Brain_ID), function(iSample){
    rowSums(exprLogic2[ ,sobj$Brain_ID == iSample])
  })
  nCells_by_donor <- do.call(cbind, nCells_by_donor)
  colnames(nCells_by_donor) <- levels(sobj$Brain_ID)
  select.genes2 <-  rownames(nCells_by_donor)[rowSums(nCells_by_donor== 0) == 0]

  sobj <- sobj[rownames(sobj) %in% select.genes2,  ]
  sobj <- sobj[rownames(sobj) %in% rowData(spe)$gene_id, ]
  # Process single cell data
  sobj <- sobj %>% 
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)
}


if(opt$data == "rat_case_control_acute" | opt$data == "rat_case_control_repeated"){
  sobj <- sce
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

  #ortholog_list_mult <- ortholog_list[!names(ortholog_list) %in% names(ortholog_list_unique)]
  #ortholog_list_mult <- lapply(ortholog_list_mult, function(iorg){
  #  if(length(iorg$rat_genes) > 1){
  #    genes2compare <- iorg$rat_genes
  #    rowmeansmat <- rowMeans(sobj[["RNA"]]$data[genes2compare,])
  #    iorg$rat_genes <- names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1]
  #  }
  #  if(length(iorg$human_genes) > 1){
  #    genes2compare <- iorg$human_genes
  #    rowmeansmat <- rowMeans(assay(spe[genes2compare, ], "logcounts"))
  #    iorg$human_genes <- names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1]
  #  }
  #  iorg
  #})
  #rat_genes <- lapply(ortholog_list_mult, function(iorg){
  #  iorg$rat_genes
  #})
  #human_genes <- lapply(ortholog_list_mult, function(iorg){
  #  iorg$human_genes
  #})
  #rat_genes <- unlist(rat_genes)
  #human_genes <- unlist(human_genes)
  #orthologs_mult_df <- data.frame("JAX_id" = names(rat_genes), "rat_genes" = rat_genes, "human_genes" = human_genes)
  #rownames(orthologs_mult_df) <- c(1:dim(orthologs_mult_df)[1])

  #orthologs_df <- rbind(orthologs_unique_df, orthologs_mult_df)
  #orthologs_df$gene_pair <- paste0(orthologs_df$rat_genes, "_", orthologs_df$human_genes)
  #orthologs_df <-  orthologs_df %>% distinct(gene_pair, .keep_all = TRUE)
  #orthologs_df <- orthologs_df[ ,c(1:3)]
  
  # Find any duplicates and only keep the one with the highest expression
  #duplicate_rat <- names(table(orthologs_df$rat_genes)[table(orthologs_df$rat_genes) > 1])
  #duplicate_human <- names(table(orthologs_df$human_genes)[table(orthologs_df$human_genes) > 1])
  #for(irat in duplicate_rat){
  #  cat(irat, "\n")
  #  df <- orthologs_df[orthologs_df$rat_genes == irat, ]
  #  genes2compare <- df$human_genes
  #  rowmeansmat <- rowMeans(assay(spe[genes2compare, ], "logcounts"))
  #  df <- df[df$human_genes == names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1], ]
  #  orthologs_df <- orthologs_df[!orthologs_df$rat_genes == irat, ]
  #  orthologs_df <- rbind(orthologs_df, df)
  #}
  #for(ihuman in duplicate_human){
  #  cat(ihuman, "\n")
  #  df <- orthologs_df[orthologs_df$human_genes == ihuman, ]
  #  genes2compare <- df$rat_genes
  #  rowmeansmat <- rowMeans(sobj[["RNA"]]$data[genes2compare,])
  #  df <- df[df$rat_genes == names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1], ]
  #  orthologs_df <- orthologs_df[!orthologs_df$human_genes == ihuman, ]
  #  orthologs_df <- rbind(orthologs_df, df)
  #}
  
  orthologs_df <- orthologs_unique_df[!duplicated(orthologs_unique_df$rat_genes), ]
  orthologs_df <- orthologs_df[!duplicated(orthologs_df$human_genes), ]
  orthologs_df <- orthologs_df[orthologs_df$rat_genes %in% rownames(sobj), ]
  orthologs_df <- orthologs_df[orthologs_df$human_genes %in% rownames(spe), ]
  saveRDS(orthologs_df, file.path(res_dir, "orthologs_df.rds"))

  # Read the differentially expressed genes in cocaine vs. saline
  diffExpr_file <- here("processed-data", "16_transfer_learning", "rat_case_control_DEGs", "NAc_rat_cocaine_all.csv")
  DEGs_cocaine <- read.csv(diffExpr_file)
  print(length(intersect(unique(DEGs_cocaine$Gene), orthologs_df$rat_genes)))

  sobj <- sobj[rownames(sobj) %in% orthologs_df$rat_genes, ]
  # Select genes for which we have orthologs
  sobj <- sobj %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)
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

  #ortholog_list_mult <- ortholog_list[!names(ortholog_list) %in% names(ortholog_list_unique)]
  #ortholog_list_mult <- lapply(ortholog_list_mult, function(iorg){
  #  if(length(iorg$rat_genes) > 1){
  #    genes2compare <- iorg$rat_genes
  #    rowmeansmat <- rowMeans(sobj[["RNA"]]$data[genes2compare,])
  #    iorg$rat_genes <- names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1]
  #  }
  #  if(length(iorg$human_genes) > 1){
  #    genes2compare <- iorg$human_genes
  #    rowmeansmat <- rowMeans(assay(spe[genes2compare, ], "logcounts"))
  #    iorg$human_genes <- names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1]
  #  }
  #  iorg
  #})
  #rat_genes <- lapply(ortholog_list_mult, function(iorg){
  #  iorg$rat_genes
  #})
  #human_genes <- lapply(ortholog_list_mult, function(iorg){
  #  iorg$human_genes
  #})
  #rat_genes <- unlist(rat_genes)
  #human_genes <- unlist(human_genes)
  #orthologs_mult_df <- data.frame("JAX_id" = names(rat_genes), "rat_genes" = rat_genes, "human_genes" = human_genes)
  #rownames(orthologs_mult_df) <- c(1:dim(orthologs_mult_df)[1])

  #orthologs_df <- rbind(orthologs_unique_df, orthologs_mult_df)
  #orthologs_df$gene_pair <- paste0(orthologs_df$rat_genes, "_", orthologs_df$human_genes)
  #orthologs_df <-  orthologs_df %>% distinct(gene_pair, .keep_all = TRUE)
  #orthologs_df <- orthologs_df[ ,c(1:3)]
  
  # Find any duplicates and only keep the one with the highest expression
  #duplicate_rat <- names(table(orthologs_df$rat_genes)[table(orthologs_df$rat_genes) > 1])
  #duplicate_human <- names(table(orthologs_df$human_genes)[table(orthologs_df$human_genes) > 1])
  #for(irat in duplicate_rat){
  #  cat(irat, "\n")
  #  df <- orthologs_df[orthologs_df$rat_genes == irat, ]
  #  genes2compare <- df$human_genes
  #  rowmeansmat <- rowMeans(assay(spe[genes2compare, ], "logcounts"))
  #  df <- df[df$human_genes == names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1], ]
  #  orthologs_df <- orthologs_df[!orthologs_df$rat_genes == irat, ]
  #  orthologs_df <- rbind(orthologs_df, df)
  #}
  #for(ihuman in duplicate_human){
  #  cat(ihuman, "\n")
  #  df <- orthologs_df[orthologs_df$human_genes == ihuman, ]
  #  genes2compare <- df$rat_genes
  #  rowmeansmat <- rowMeans(sobj[["RNA"]]$data[genes2compare,])
  #  df <- df[df$rat_genes == names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1], ]
  #  orthologs_df <- orthologs_df[!orthologs_df$human_genes == ihuman, ]
  #  orthologs_df <- rbind(orthologs_df, df)
  #}
  
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


plot_dir <- paste0(plot_dir, "/", opt$data)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
if(opt$data == "human_NAc"){
pdf(file.path(plot_dir, "snRNA_seq_PCs.pdf"), height = 4, width = 6)
# Check Elbow plot to determine relevant number of PCs
ElbowPlot(sobj, ndims = 50)
DimPlot(sobj, reduction = "pca", group.by = "Sample")
DimPlot(sobj, reduction = "pca", group.by = "Brain_ID")
DimPlot(sobj, reduction = "pca", group.by = "Sort")
DimPlot(sobj, reduction = "pca", group.by = "snRNA_date")
FeaturePlot(sobj, reduction = "pca", features = c("nCount_RNA"))
FeaturePlot(sobj, reduction = "pca", features = c("nFeature_RNA"))
FeaturePlot(sobj, reduction = "pca", features = c("percent_mito"))
dev.off()
}
if(opt$data == "rat_case_control_acute" | opt$data == "rat_case_control_repeated"){
pdf(file.path(plot_dir, "snRNA_seq_PCs.pdf"), height = 4, width = 6)
# Check Elbow plot to determine relevant number of PCs
ElbowPlot(sobj, ndims = 50)
DimPlot(sobj, reduction = "pca", group.by = "Sex")
DimPlot(sobj, reduction = "pca", group.by = "GEM")
DimPlot(sobj, reduction = "pca", group.by = "Stim")
DimPlot(sobj, reduction = "pca", group.by = "CellType")
DimPlot(sobj, reduction = "pca", group.by = "Combo_CellType")
FeaturePlot(sobj, reduction = "pca", features = c("nCount_RNA"))
FeaturePlot(sobj, reduction = "pca", features = c("nFeature_RNA"))
FeaturePlot(sobj, reduction = "pca", features = c("percent_mito"))
dev.off()
}
if(opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
  pdf(file.path(plot_dir, "snRNA_seq_PCs.pdf"), height = 4, width = 6)
  # Check Elbow plot to determine relevant number of PCs
  ElbowPlot(sobj, ndims = 50)
  DimPlot(sobj, reduction = "pca", group.by = "subject")
  DimPlot(sobj, reduction = "pca", group.by = "treatment")
  DimPlot(sobj, reduction = "pca", group.by = "CellType")
  FeaturePlot(sobj, reduction = "pca", features = c("nCount_RNA"))
  FeaturePlot(sobj, reduction = "pca", features = c("nFeature_RNA"))
  FeaturePlot(sobj, reduction = "pca", features = c("percent_mito"))
  dev.off()
}


# Get the proportion of variance of each PC explained by meta-data
pca_embeddings <- Embeddings(sobj, reduction = "pca")[, 1:20]

if(opt$data == "human_NAc"){
  pca_res <- matrix(NA, nrow = dim(pca_embeddings)[2], ncol = 11)
for(i in c(1:dim(pca_embeddings)[2])){
  cat(i, "\n")
  pca_res[i, 1] <- summary(lm(pca_embeddings[ ,i] ~ sobj$nCount_RNA))$adj.r.squared
  pca_res[i, 2] <- summary(lm(pca_embeddings[ ,i] ~ sobj$nFeature_RNA))$adj.r.squared
  pca_res[i, 3] <- summary(lm(pca_embeddings[ ,i] ~ sobj$percent_mito))$adj.r.squared
  pca_res[i, 4] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$Sample)))$adj.r.squared
  pca_res[i, 5] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$Brain_ID)))$adj.r.squared
  pca_res[i, 6] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$CellType)))$adj.r.squared
  pca_res[i, 7] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$Sort)))$adj.r.squared
  pca_res[i, 8] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$snRNA_date)))$adj.r.squared
  pca_res[i, 9] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$Chromium_cDNA_date)))$adj.r.squared
  pca_res[i, 10] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$Chromium_library_date)))$adj.r.squared
  pca_res[i, 11] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$Chromium_kit)))$adj.r.squared
}
colnames(pca_res) <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "Sample", "Brain_ID", "Cell_Type", "Sort", "snRNA_date", "Chromium_cDNA_date", "Chromium_library_date", "Chromium_kit")
rownames(pca_res) <- paste0("PC_", c(1:20))
pca_res <- reshape2::melt(pca_res)
}

if(opt$data == "rat_case_control_acute" | opt$data == "rat_case_control_repeated"){
   pca_res <- matrix(NA, nrow = dim(pca_embeddings)[2], ncol = 8)
for(i in c(1:dim(pca_embeddings)[2])){
  cat(i, "\n")
  pca_res[i, 1] <- summary(lm(pca_embeddings[ ,i] ~ sobj$nCount_RNA))$adj.r.squared
  pca_res[i, 2] <- summary(lm(pca_embeddings[ ,i] ~ sobj$nFeature_RNA))$adj.r.squared
  pca_res[i, 3] <- summary(lm(pca_embeddings[ ,i] ~ sobj$percent_mito))$adj.r.squared
  pca_res[i, 4] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$Stim)))$adj.r.squared
  pca_res[i, 5] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$Sex)))$adj.r.squared
  pca_res[i, 6] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$GEM)))$adj.r.squared
  pca_res[i, 7] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$CellType)))$adj.r.squared
  pca_res[i, 8] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$Combo_CellType)))$adj.r.squared
}
colnames(pca_res) <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "Stim", "Sex", "GEM", "CellType", "Combo_CellType")
rownames(pca_res) <- paste0("PC_", c(1:20))
pca_res <- reshape2::melt(pca_res)
}

if(opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
   pca_res <- matrix(NA, nrow = dim(pca_embeddings)[2], ncol = 6)
for(i in c(1:dim(pca_embeddings)[2])){
  cat(i, "\n")
  pca_res[i, 1] <- summary(lm(pca_embeddings[ ,i] ~ sobj$nCount_RNA))$adj.r.squared
  pca_res[i, 2] <- summary(lm(pca_embeddings[ ,i] ~ sobj$nFeature_RNA))$adj.r.squared
  pca_res[i, 3] <- summary(lm(pca_embeddings[ ,i] ~ sobj$percent_mito))$adj.r.squared
  pca_res[i, 4] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$subject)))$adj.r.squared
  pca_res[i, 5] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$treatment)))$adj.r.squared
  pca_res[i, 6] <- summary(lm(pca_embeddings[ ,i] ~ as.factor(sobj$CellType)))$adj.r.squared
}
colnames(pca_res) <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "Subject", "Treatment", "CellType")
rownames(pca_res) <- paste0("PC_", c(1:20))
pca_res <- reshape2::melt(pca_res)
}

pdf(file.path(plot_dir, "snRNA_seq_PC_drivers.pdf"), height = 6, width = 8)
ggplot(pca_res, aes(x = Var1, y = Var2, fill = value)) + geom_tile(color = "black") + 
scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) + coord_fixed() + xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()

sobj <- sobj %>%
  FindNeighbors(dims = 1:20, verbose = FALSE) %>%
  FindClusters(resolution = 1.0, verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)

saveRDS(sobj, file.path(res_dir, "snRNA_seq_NAc.rds"))

if(opt$data == "human_NAc"){
  pdf(file.path(plot_dir, "UMAP.pdf"), height = 6, width = 8)
  DimPlot(sobj, group.by = "CellType")
  dev.off()
}

if(opt$data == "rat_case_control_acute" | opt$data == "rat_case_control_repeated"){
  pdf(file.path(plot_dir, "UMAP.pdf"), height = 6, width = 8)
  DimPlot(sobj, group.by = "Combo_CellType")
  DimPlot(sobj, group.by = "Stim")
  dev.off()
}

if(opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
  pdf(file.path(plot_dir, "UMAP.pdf"), height = 6, width = 8)
  DimPlot(sobj, group.by = "CellType")
  dev.off()
}

sessionInfo()