library(here)
library(SingleCellExperiment)
library(jaffelab)
library(scater)
library(scran)
library(readxl)
library(Polychrome)
library(cluster)
library(limma)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(ggrastr)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(spatialNAcUtils)
library(Polychrome)
data(palette36)
library(getopt)

# Read in settings
spec <- matrix(
    c("capture_area_path", "c", 1, "character", "Specify the location of capture area results", 
      "donor_path", "d", 1, "character", "Specify the location of donor results", 
       "cluster_col", "p", 1, "character", "Specify cluster column", 
       "spot_path", "s", 1, "character", "Specify location of spot results"),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)
#opt <- list()
#opt$capture_area_path <- '02_BayesSpace/pseudobulk_capture_area'
#opt$donor_path <- '01_precast/pseudobulk_donor/final_clusters'
#opt$cluster_col <- 'BayesSpace_harmony_k3'
#opt$spot_path <- '01_precast/nnSVG_precast/random_start_1/PRECAST_k3_marker_genes.csv'

if(opt$capture_area_path == "NA"){
  opt$capture_area_path <- NA
}

if(opt$donor_path == "NA"){
  opt$donor_path <- NA
}

if(opt$spot_path == "NA"){
  opt$spot_path <- NA
}


if(grepl("BayesSpace", opt$cluster_col)){
    clust_col_n <- as.numeric(gsub("BayesSpace_harmony_k", "", opt$cluster_col))
    clust_col_n_nice <- sprintf("%02d", clust_col_n)
    opt$cluster_col <- paste0("BayesSpace_harmony_k", clust_col_n_nice)
}
# Plot the enhanced volcano plots
if(!is.na(opt$spot_path)){
  spot_level_marker_fn <- here('processed-data', '07_spatial_domains', opt$spot_path)
}
if(!is.na(opt$donor_path)){
  pseudobulk_donor_marker_fn <- here('processed-data', '10_post_clustering_analysis', '01_pseudobulk_markers', 
  opt$donor_path ,paste0("model_results_",opt$cluster_col, "_FDR5perc.csv"))
}
if(!is.na(opt$capture_area_path)){
  pseudobulk_capture_area_marker_fn <- here('processed-data', '10_post_clustering_analysis', '01_pseudobulk_markers',
  opt$capture_area_path, paste0("model_results_",opt$cluster_col ,"_FDR5perc.csv"))
}

if(grepl("precast", opt$capture_area_path)){
  plotDir <- here('plots', '10_post_clustering_analysis', '01_pseudobulk_markers', unlist(strsplit(opt$capture_area_path, split = "/"))[1],
 'visualize_markers', unlist(strsplit(opt$capture_area_path, split = "/"))[3], opt$cluster_col)
}

if(grepl("BayesSpace", opt$capture_area_path)){
  plotDir <- here('plots', '10_post_clustering_analysis', '01_pseudobulk_markers', unlist(strsplit(opt$capture_area_path, split = "/"))[1],
  'visualize_markers' , opt$cluster_col)
}

dir.create(plotDir, recursive = TRUE)

# Read in the enrichment results and make plots
# Spot level
if(!is.na(opt$spot_path)){
  spot_level_markers <- read.csv(spot_level_marker_fn)
  pVolcano_spot <- list()
  for(clust in unique(spot_level_markers$cluster)){
    spot_enr <- spot_level_markers[spot_level_markers$cluster == clust, ]
    pVolcano_spot[[clust]] <- EnhancedVolcano(spot_enr, lab = spot_enr$gene, x = 'avg_log2FC', y = 'p_val_adj', title = paste0("Cluster_", clust),
    pointSize = 1.0,labSize = 3.0, colAlpha = 0.5, pCutoff = 1e-3, FCcutoff = 1, legendLabels=c('Not sig.','Log FC','p-value',
      'p-value & Log FC'), drawConnectors = TRUE, min.segment.length = 0.5) 
  }
  pdf(file.path(plotDir, "enhanced_volcano_spot.pdf"), width = 10, height = 10)
  for(k in 1:length(pVolcano_spot)){
    print(pVolcano_spot[[k]])
    }
  dev.off() 
}

# Donor level
if(!is.na(opt$donor_path)){
  pseudobulk_donor_markers <- read.csv(pseudobulk_donor_marker_fn)
  donor_enr <- pseudobulk_donor_markers[pseudobulk_donor_markers$model_type == "enrichment", ]
  pVolcano_donor <- list()
  for(clust in unique(donor_enr$test)){
    print(clust)
    donor_enr_clust <- donor_enr[donor_enr$test == clust, ]
    pVolcano_donor[[clust]] <- EnhancedVolcano(donor_enr_clust, lab = donor_enr_clust$gene, x = 'logFC', y = 'fdr', title = clust,
    pointSize = 1.0,labSize = 3.0, colAlpha = 0.5, pCutoff = 1e-3, FCcutoff = 1, legendLabels=c('Not sig.','Log FC','p-value',
      'p-value & Log FC'), drawConnectors = TRUE, min.segment.length = 0.5, max.overlaps = 5) 
  }
  pdf(file.path(plotDir, "enhanced_volcano_donor.pdf"), width = 10, height = 10)
  for(k in 1:length(pVolcano_donor)){
    print(pVolcano_donor[[k]])
  }
  dev.off() 
}

# Capture area level
if(!is.na(opt$capture_area_path)){
  pseudobulk_capture_area_markers <- read.csv(pseudobulk_capture_area_marker_fn)
  capture_area_enr <- pseudobulk_capture_area_markers[pseudobulk_capture_area_markers$model_type == "enrichment", ]
  pVolcano_captureArea <- list()
  for(clust in unique(capture_area_enr$test)){
    print(clust)
    capture_area_enr_clust <- capture_area_enr[capture_area_enr$test == clust, ]
    pVolcano_captureArea[[clust]] <- EnhancedVolcano(capture_area_enr_clust, lab = capture_area_enr_clust$gene, x = 'logFC', y = 'fdr', title = clust,
    pointSize = 1.0,labSize = 3.0, colAlpha = 0.5, pCutoff = 1e-3, FCcutoff = 1, legendLabels=c('Not sig.','Log FC','p-value',
      'p-value & Log FC'), drawConnectors = TRUE, min.segment.length = 0.5, max.overlaps = 5) 
  }
  pdf(file.path(plotDir, "enhanced_volcano_captureArea.pdf"), width = 10, height = 10)
  for(k in 1:length(pVolcano_captureArea)){
    print(pVolcano_captureArea[[k]])
  }
  dev.off() 
}

# Compare across the different methods used in calling DE genes if all three methods are present
if(!is.na(opt$capture_area_path) & !is.na(opt$donor_path) & !is.na(opt$spot_path)){
  pHits <- list()
  pOverlap <- list()
  for(clust in unique(spot_level_markers$cluster)){
    spot_enr <- spot_level_markers[spot_level_markers$cluster == clust, ]
    donor_enr <- pseudobulk_donor_markers[pseudobulk_donor_markers$model_type == "enrichment", ]
    donor_enr <- donor_enr[donor_enr$test == paste0("X", clust), ]
    capture_area_enr <- pseudobulk_capture_area_markers[pseudobulk_capture_area_markers$model_type == "enrichment", ]
    capture_area_enr <- capture_area_enr[capture_area_enr$test == paste0("X", clust), ]

    spot_enr <- spot_enr[spot_enr$p_val_adj < 1e-3, ]
    capture_area_enr <- capture_area_enr[capture_area_enr$fdr < 1e-3, ]
    donor_enr <- donor_enr[donor_enr$fdr < 1e-3, ]

    spot_enr <- spot_enr[order(-spot_enr$avg_log2FC), ]
    donor_enr <- donor_enr[order(-donor_enr$logFC), ]
    capture_area_enr <- capture_area_enr[order(-capture_area_enr$logFC), ]

    min_nHits <- min(dim(spot_enr)[1], dim(donor_enr)[1], dim(capture_area_enr)[1])
    if(min_nHits == 0){
      next
    }
    nHits <- data.frame("Method" = c("Spot-level", "Capture-area level", "Donor-level"), 
                        "nHits_FDR_0.001" = c(dim(spot_enr)[1], dim(capture_area_enr)[1], dim(donor_enr)[1]))
    nHits$Method <- factor(nHits$Method, levels = c("Spot-level", "Capture-area level", "Donor-level"))
  
    pHits[[clust]] <- ggplot(nHits, aes(x = Method, y = nHits_FDR_0.001, fill = Method)) + geom_bar(stat="identity") + ggtitle(paste0("Cluster_", clust))

    T <- min(dim(spot_enr)[1], dim(donor_enr)[1], dim(capture_area_enr)[1])
    ind <- 1
    spot_and_donor <- c()
    donor_and_capture_area <- c()
    spot_and_capture_area <- c()
    all_methods <- c()
    if(T < 100){
      stepSize <- 20
    }
    if(T >= 100 & T < 500){
      stepSize <- 10
    }
    if(T >= 500){
      stepSize <- 8
    }
    for(t in as.integer(seq(T/stepSize, T, T/stepSize))){
        spot_and_donor[ind] <- length(intersect(spot_enr$gene[1:t], donor_enr$gene[1:t]))
        donor_and_capture_area[ind] <- length(intersect(donor_enr$gene[1:t], capture_area_enr$gene[1:t]))
        spot_and_capture_area[ind] <- length(intersect(spot_enr$gene[1:t], capture_area_enr$gene[1:t]))
        all_methods[ind] <- length(intersect(spot_enr$gene[1:t], intersect(spot_enr$gene[1:t], capture_area_enr$gene[1:t])))
        ind <- ind+1
    }

    summary_res <- data.frame("Cutoffs" = as.integer(seq(T/stepSize, T, T/stepSize)), "Spot_and_donor" = spot_and_donor, "Donor_and_capture_area" = donor_and_capture_area, 
    "Spot_and_capture_area" = spot_and_capture_area, "All_methods" = all_methods)
    summary_res$Spot_and_donor <- summary_res$Spot_and_donor/summary_res$Cutoffs
    summary_res$Donor_and_capture_area <- summary_res$Donor_and_capture_area/summary_res$Cutoffs
    summary_res$Spot_and_capture_area <- summary_res$Spot_and_capture_area/summary_res$Cutoffs
    summary_res$All_methods <- summary_res$All_methods/summary_res$Cutoffs
    summary_res <- reshape2::melt(summary_res, id.vars = "Cutoffs")
    colnames(summary_res) <- c("Cutoffs", "Comparison", "Fraction_overlap")
    summary_res$Comparison <- factor(summary_res$Comparison, levels = c("Donor_and_capture_area", "Spot_and_donor", "Spot_and_capture_area", "All_methods"))

    pOverlap[[clust]] <- ggplot(summary_res, aes(x = Cutoffs, y = Fraction_overlap, color = Comparison, fill = Comparison)) + geom_bar(stat = "identity", position = "dodge") + ggtitle(paste0("Cluster_", clust))
  }
  pdf(file.path(plotDir, "overlap_across_methods.pdf"), width = 8, height = 4)
  for(k in 1:length(pOverlap)){
    print(pOverlap[[k]] + theme_classic())
  }
  dev.off()

  pdf(file.path(plotDir, "nHits_across_methods.pdf"), width = 6, height = 4)
  for(k in 1:length(pHits)){
    print(pHits[[k]]+ theme_classic())
  }
  dev.off()
}
