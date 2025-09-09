rm(list = ls())
library(Seurat)
library(Matrix)
library(CoGAPS)
library(getopt)
library(here)
library(ggplot2)
library(pheatmap)
library(HDF5Array)
library(SpatialExperiment)
library(spatialLIBD)
library(scran)
library(scuttle)
library(SingleCellExperiment)
library(RcppML)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggsignif)

# Define command-line arguments
spec <- matrix(c(
  'data',        'd', 1, 'character', 'Specify the input dataset',
  'nFactors',    'n', 1, 'integer',   'Number of NMF factors to use',
  'select_HVGs', 's', 1, 'logical',   'Whether to use highly variable genes (TRUE/FALSE)'
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

# Check if all required arguments were provided
if (is.null(opt$data) || is.null(opt$nFactors) || is.null(opt$select_HVGs)) {
  cat(getopt(spec, usage = TRUE))
  stop("Missing required arguments.\n")
}

# Display the selected options
cat("Selected options:\n")
print(opt)

#opt <- list()
#opt$data <- "rat_case_control_morphine_acute"
#opt$nFactors <- 30
#opt$select_HVGs <- FALSE

# Set the locations of the output files and plots
res_dir <- here::here("processed-data", "16_transfer_learning", "04_cellType_NMF")
plot_dir <- here::here("plots", "16_transfer_learning", "04_cellType_NMF")
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
select.genes <-  readRDS(file.path(here::here("processed-data", "16_transfer_learning", "04_cellType_NMF"), "select_genes.rds"))
# Only select those genes which have non-zero expression in atleast 1 spot in each slide
spe <- spe[rownames(spe) %in% select.genes, ]
spe <- logNormCounts(spe)

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

# Subset only to the MSN cell types in the data
if(opt$data == "rat_case_control_cocaine_acute" | opt$data == "rat_case_control_cocaine_repeated"){
  sce$Combo_CellType <- as.character(sce$Combo_CellType)
  sce$is_MSN <- sce$Combo_CellType %in% c("Drd1-MSN-1", "Drd1-MSN-2", "Drd2-MSN-1", "Drd2-MSN-2", "Drd3-MSN", "Grm8-MSN")
  sce <- subset(sce, is_MSN == TRUE)
  nonzero_genes <- rowSums(sce@assays$RNA@counts) > 0
  # Subset the Seurat object to keep only genes with non-zero expression
  sce <- subset(sce, features = names(nonzero_genes[nonzero_genes]))  
  

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

  # Read the differentially expressed genes in cocaine vs. saline
  diffExpr_file <- here("processed-data", "16_transfer_learning", "rat_case_control_DEGs", "NAc_rat_cocaine_all.csv")
  DEGs_cocaine <- read.csv(diffExpr_file)
  print(length(intersect(unique(DEGs_cocaine$Gene), orthologs_df$rat_genes)))

  sce <- sce[rownames(sce) %in% orthologs_df$rat_genes, ]
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  # Select genes for which we have orthologs
  sce <- sce %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)

  sce <- sce %>%
  FindNeighbors(dims = 1:20, verbose = FALSE) %>%
  FindClusters(resolution = 1.0, verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)
  
  pdf(file.path(plot_dir, "UMAP.pdf"))
  DimPlot(sce, group.by = "Combo_CellType")
  DimPlot(sce, group.by = "Stim")
  dev.off() 
  
  if(opt$select_HVGs){
    hvg_genes <- VariableFeatures(sce)
    sce2 <- subset(sce, features = hvg_genes)
    sce2 <- NormalizeData(sce2)
    dat <- sce2[["RNA"]]$data
    dat <- dat[rowSums(dat) > 0, ]
  }else{
    dat <- sce[["RNA"]]$data
    dat <- dat[rowSums(dat) > 0, ]
  }

  nmf_precomputed <- TRUE
  if(nmf_precomputed){
    if(opt$select_HVGs){
      nmf_results <- readRDS(file.path(res_dir, paste0("nmf_results_HVGs_", opt$nFactors, "factors.rds")))
    }else{
      nmf_results <- readRDS(file.path(res_dir, paste0("nmf_results_", opt$nFactors, "factors.rds")))
    }
  }else{
      # Run NMF
    nmf_results <- RcppML::nmf(dat,
                 k=opt$nFactors,
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 1135,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE)
    if(opt$select_HVGs){
        saveRDS(nmf_results, file.path(res_dir, paste0("nmf_results_HVGs_", opt$nFactors, "factors.rds")))
    }else{
        saveRDS(nmf_results, file.path(res_dir, paste0("nmf_results_", opt$nFactors, "factors.rds")))
    }
  }
  
  sce@meta.data <- cbind(sce@meta.data, data.frame(t(nmf_results@h)))
  sce$CellType_Stim <- paste0(sce$Combo_CellType, "_", sce$Stim)
  df <- sce@meta.data
  nmf_cols <- grep("^nmf", colnames(df), value = TRUE)
  # Step 2: Reshape to long format
  nmf_long <- df %>%
    select(Combo_CellType, Stim, CellType_Stim, all_of(nmf_cols)) %>%
    pivot_longer(cols = all_of(nmf_cols),
                names_to = "NMF_Component",
                values_to = "Weight")

  # Step 3: Run Wilcoxon tests
  effect_results <- nmf_long %>%
    group_by(Combo_CellType, NMF_Component) %>%
    wilcox_test(Weight ~ Stim, detailed = TRUE) %>%
    mutate(Adj_P = p.adjust(p, method = "BH")) %>%
    add_significance("Adj_P") %>%
    rename(sig = Adj_P.signif)

  # Step 4: Compute effect sizes (fast version)
  effect_sizes <- nmf_long %>%
    group_by(Combo_CellType, NMF_Component) %>%
    wilcox_effsize(Weight ~ Stim)

  # Step 5: Merge test + effect size results
  effect_results <- left_join(effect_results, effect_sizes,
                              by = c("Combo_CellType", "NMF_Component"))

  # Step 6: Filter significant and biologically meaningful effects
  filtered_effects <- effect_results %>%
    filter(Adj_P < 0.05, effsize > 0.3)

  # Step 7: Merge with long data for plotting
  plot_data <- nmf_long %>%
    left_join(effect_results %>% select(Combo_CellType, NMF_Component, Adj_P, sig), 
              by = c("Combo_CellType", "NMF_Component"))
  
  # Determine plotting order: significant first, then non-significant
  all_components <- unique(nmf_long$NMF_Component)
  significant_components <- unique(filtered_effects$NMF_Component)
  non_significant_components <- setdiff(all_components, significant_components)
  ordered_components <- c(significant_components, non_significant_components)

  pdf(file.path(plot_dir, paste0("nmf_VlnPlot_", 
    ifelse(opt$select_HVGs, "HVGs_", ""), 
    opt$nFactors, "_factors.pdf")), width = 10, height = 5)

  # Loop through all components in desired order
  for (nmf in ordered_components) {
    nmf_data <- plot_data %>% filter(NMF_Component == nmf)
    nmf_sig <- filtered_effects %>% filter(NMF_Component == nmf)
  
    # Default plot
    p <- ggplot(nmf_data, aes(x = CellType_Stim, y = Weight, fill = Stim)) +
      geom_violin(scale = "width", trim = FALSE) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.2) +
      labs(title = paste("NMF Component:", nmf),
          subtitle = ifelse(nrow(nmf_sig) > 0, "Significance: *, **, *** with effect size", "Not significant"),
          x = "Cell Type + Treatment (CellType_Stim)",
          y = "NMF Weight") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold"))

    # Add significance annotations if applicable
    if (nrow(nmf_sig) > 0) {
      nmf_sig <- nmf_sig %>%
        mutate(group1 = paste(Combo_CellType, "Cocaine", sep = "_"),
              group2 = paste(Combo_CellType, "Saline", sep = "_"),
              label = paste0(sig, " (", round(effsize, 2), ")"))

      comparison_pairs <- split(nmf_sig[, c("group1", "group2")], seq(nrow(nmf_sig)))
      comparison_pairs <- lapply(comparison_pairs, function(row) c(row$group1, row$group2))

      y_pos_map <- nmf_sig %>%
        rowwise() %>%
        mutate(y_pos = max(nmf_data$Weight[nmf_data$CellType_Stim %in% c(group1, group2)], na.rm = TRUE) * 1.5) %>%
        pull(y_pos)

      p <- p + geom_signif(
        comparisons = comparison_pairs,
        annotations = nmf_sig$label,
        y_position = y_pos_map,
        tip_length = 0.01,
        textsize = 3,
        vjust = 0
      )
    }
    print(p)
  }
  dev.off()
}

if(opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
  sce$cluster <- as.character(sce$cluster)
  sce$is_MSN <- sce$cluster %in% c("Drd1 MSN1", "Drd1 MSN2", "Drd1 MSN3","Drd2 MSN1", "Drd2 MSN2", "Drd2 MSN3", 
  "Drd2 MSN4", "Grm8 MSN", "MSN1", "MSN2", "MSN3")
  sce <- sce[ ,sce$is_MSN]
  nonzero_genes <- rowSums(counts(sce)) > 0
  sce <- sce[nonzero_genes, ]
  rownames(sce) <- make.unique(rowData(sce)$gene_name)
  
  # Create a seurat object
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
  sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)
  # Select genes for which we have orthologs
  sobj <- sobj %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)

  sobj <- sobj %>%
  FindNeighbors(dims = 1:20, verbose = FALSE) %>%
  FindClusters(resolution = 1.0, verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)
  
  pdf(file.path(plot_dir, "UMAP.pdf"))
  DimPlot(sobj, group.by = "CellType")
  DimPlot(sobj, group.by = "treatment")
  dev.off() 

  if(opt$select_HVGs){
    hvg_genes <- VariableFeatures(sobj)
    sobj2 <- subset(sobj, features = hvg_genes)
    sobj2 <- NormalizeData(sobj2)
    dat <- sobj2[["RNA"]]$data
    dat <- dat[rowSums(dat) > 0, ]
  }else{
    dat <- sobj[["RNA"]]$data
    dat <- dat[rowSums(dat) > 0, ]
  }
  nmf_precomputed <- TRUE
  if(nmf_precomputed){
    if(opt$select_HVGs){
        nmf_results <- readRDS(file.path(res_dir, paste0("nmf_results_HVGs_", opt$nFactors, "factors.rds")))
    }else{
        nmf_results <- readRDS(file.path(res_dir, paste0("nmf_results_", opt$nFactors, "factors.rds")))
    }
  }else{
      # Run NMF
    nmf_results <- RcppML::nmf(dat,
                 k=opt$nFactors,
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 1135,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE)
    if(opt$select_HVGs){
        saveRDS(nmf_results, file.path(res_dir, paste0("nmf_results_HVGs_", opt$nFactors, "factors.rds")))
    }else{
        saveRDS(nmf_results, file.path(res_dir, paste0("nmf_results_", opt$nFactors, "factors.rds")))
    }

  }
  
  sobj@meta.data <- cbind(sobj@meta.data, data.frame(t(nmf_results@h)))
  sobj$CellType_Stim <- paste0(sobj$CellType, "_", sobj$treatment)
  df <- sobj@meta.data
  nmf_cols <- grep("^nmf", colnames(df), value = TRUE)
  # Step 2: Reshape to long format
  nmf_long <- df %>%
    select(CellType, treatment, CellType_Stim, all_of(nmf_cols)) %>%
    pivot_longer(cols = all_of(nmf_cols),
                names_to = "NMF_Component",
                values_to = "Weight")

  # Step 3: Run Wilcoxon tests
  effect_results <- nmf_long %>%
    group_by(CellType, NMF_Component) %>%
    wilcox_test(Weight ~ treatment, detailed = TRUE) %>%
    mutate(Adj_P = p.adjust(p, method = "BH")) %>%
    add_significance("Adj_P") %>%
    dplyr::rename(sig = Adj_P.signif)

  # Step 4: Compute effect sizes (fast version)
  effect_sizes <- nmf_long %>%
    group_by(CellType, NMF_Component) %>%
    wilcox_effsize(Weight ~ treatment)

  # Step 5: Merge test + effect size results
  effect_results <- left_join(effect_results, effect_sizes,
                              by = c("CellType", "NMF_Component"))

  # Step 6: Filter significant and biologically meaningful effects
  filtered_effects <- effect_results %>%
    filter(Adj_P < 0.05, effsize > 0.3)

  # Step 7: Merge with long data for plotting
  plot_data <- nmf_long %>%
    left_join(effect_results %>% select(CellType, NMF_Component, Adj_P, sig), 
              by = c("CellType", "NMF_Component"))

  # Determine plotting order: significant first, then non-significant
  all_components <- unique(nmf_long$NMF_Component)
  significant_components <- unique(filtered_effects$NMF_Component)
  non_significant_components <- setdiff(all_components, significant_components)
  ordered_components <- c(significant_components, non_significant_components)

  pdf(file.path(plot_dir, paste0("nmf_VlnPlot_", 
    ifelse(opt$select_HVGs, "HVGs_", ""), 
    opt$nFactors, "_factors_1.pdf")), width = 10, height = 5)
  # Loop through all components in desired order
  for (nmf in ordered_components) {
    nmf_data <- plot_data %>% filter(NMF_Component == nmf)
    nmf_sig <- filtered_effects %>% filter(NMF_Component == nmf)
  
    # Default plot
    p <- ggplot(nmf_data, aes(x = CellType_Stim, y = Weight, fill = treatment)) +
      geom_violin(scale = "width", trim = FALSE) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.2) +
      labs(title = paste("NMF Component:", nmf),
          subtitle = ifelse(nrow(nmf_sig) > 0, "Significance: *, **, *** with effect size", "Not significant"),
          x = "Cell Type + Treatment (CellType_Stim)",
          y = "NMF Weight") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold"))

    # Add significance annotations if applicable
    if (nrow(nmf_sig) > 0) {
      nmf_sig <- nmf_sig %>%
        mutate(group1 = paste(CellType, "morphine", sep = "_"),
              group2 = paste(CellType, "saline", sep = "_"),
              label = paste0(sig, " (", round(effsize, 2), ")"))

      comparison_pairs <- split(nmf_sig[, c("group1", "group2")], seq(nrow(nmf_sig)))
      comparison_pairs <- lapply(comparison_pairs, function(row) c(row$group1, row$group2))

      y_pos_map <- nmf_sig %>%
        rowwise() %>%
        mutate(y_pos = max(nmf_data$Weight[nmf_data$CellType_Stim %in% c(group1, group2)], na.rm = TRUE) * 1.5) %>%
        pull(y_pos)

      p <- p + geom_signif(
        comparisons = comparison_pairs,
        annotations = nmf_sig$label,
        y_position = y_pos_map,
        tip_length = 0.01,
        textsize = 3,
        vjust = 0
      )
    }
    print(p)
  }
  dev.off()
}
