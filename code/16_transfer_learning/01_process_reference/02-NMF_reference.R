library(RcppML)
library(pheatmap)
library(SingleCellExperiment)
library(here)
library(scuttle)
library(dplyr)
library(getopt)
library(Matrix)
library(HDF5Array)
set.seed(1234)

spec <- matrix(
    c(
        "gene_selection_strategy", "g", 1, "character", "Choose all genes, or highly deviant genes based on snRNA-seq, or nnSVGs", 
        "data", "d", 1, "character", "Specify the input dataset"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

#opt <- list()
#opt$gene_selection_strategy <- "all_genes"
#opt$data <- "human_NAc"
print(opt$gene_selection_strategy)
print(opt$data)

# Read data based on whether we are processing the human or rat NAc snRNA seq data
dat_dir <- here::here("processed-data", "12_snRNA")
res_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "RCppML", opt$data)
plot_dir <- here::here("plots", "16_transfer_learning", "01_process_reference", "RCppML", opt$data)
dir.create(res_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

if(opt$data == "human_NAc"){
  sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
}else{
  if(opt$data == "rat_case_control"){
    sce <- readRDS(file = file.path(dat_dir, "NAc_Combo_Integrated.RDS"))
  }else{
    stop("Invalid input data set")
  }
}

if(opt$data == "human_NAc"){
    if(opt$gene_selection_strategy == "all_genes"){
      # Read in the spatial experiment data
      raw_in_path <- here("processed-data", "05_harmony_BayesSpace", "02-compute_QC_metrics", "spe_with_QC_metrics_hdf5")
      spe <- loadHDF5SummarizedExperiment(raw_in_path)
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
      common_genes <- intersect(rownames(spe), rownames(sce))
      sce <- sce[rownames(sce) %in% common_genes, ]
    }else{
      if(opt$gene_selection_strategy == "snRNA_seq"){
      select_genes <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing = TRUE)][1:5000]
      sce <- sce[ ,rownames(sce) %in% select_genes]
      }else{
        if(opt$gene_selection_strategy == "nnSVG"){
          svg_path <- here(
          "processed-data", "05_harmony_BayesSpace", "07-run_nnSVG", "nnSVG_precast_out",
          "summary_across_samples.csv")
          num_genes <- 5000
          select_genes <- read.csv(svg_path) |>
            as_tibble() |>
            arrange(nnsvg_avg_rank_rank) |>
            slice_head(n = num_genes) |>
            pull(gene_id)
            sce <- sce[ ,rownames(sce) %in% select_genes]
        }else{
          stop("Invalid gene selection strategy")
      }
  }
}
}else{
  if(opt$data == "rat_case_control"){
    if(opt$gene_selection_strategy == "all_genes"){
      refDir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "preliminary_analysis")
      orthologs_df <- readRDS(file.path(refDir, opt$data, "orthologs_df.rds"))
      sce <- sce[rownames(sce) %in% orthologs_df$rat_genes, ]
    }else{
      if(opt$gene_selection_strategy == "snRNA_seq"){
        refDir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "preliminary_analysis")
        orthologs_df <- readRDS(file.path(refDir, opt$data, "orthologs_df.rds"))
        sce <- sce[rownames(sce) %in% orthologs_df$rat_genes, ]
        sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000) 
        var_genes <- VariableFeatures(sce)
        select_genes <- intersect(orthologs_df$rat_genes, var_genes)
        sce <- sce[rownames(sce) %in% select_genes, ]
      }else{
        if(opt$gene_selection_strategy == "nnSVG"){
          refDir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "preliminary_analysis")
          orthologs_df <- readRDS(file.path(refDir, opt$data, "orthologs_df.rds"))
          svg_path <- here(
          "processed-data", "05_harmony_BayesSpace", "07-run_nnSVG", "nnSVG_precast_out",
          "summary_across_samples.csv")
          select_genes <- read.csv(svg_path)|>
            arrange(nnsvg_avg_rank_rank)
          raw_in_path <- here("processed-data", "05_harmony_BayesSpace", "02-compute_QC_metrics", "spe_with_QC_metrics_hdf5")
          spe <- loadHDF5SummarizedExperiment(raw_in_path)
          geneData <- rowData(spe)
          geneData <- geneData[rownames(geneData) %in% select_genes$gene_id, ]
          geneData <- geneData[match(select_genes$gene_id, rownames(geneData)), ]
          select_genes$gene_name <- geneData$gene_name
          select_genes <- select_genes[select_genes$gene_name %in% orthologs_df$human_genes, ]
          num_genes <- 5000
          select_genes <- select_genes[1:min(num_genes, dim(select_genes)[1]), "gene_name"]
          orthologs_df <- orthologs_df[orthologs_df$human_genes %in% select_genes, ]
          rat_select_genes <- orthologs_df$rat_genes
          sce <- sce[rownames(sce) %in% rat_select_genes, ]
        }else{
          stop("Invalid gene selection")
        }
      }
    }
  }
}

if(opt$data == "human_NAc"){
  dat <- assay(sce,'logcounts')
}
if(opt$data == "rat_case_control"){
  dat <- sce[["RNA"]]$data
}


print("Running RCppML")
options(RcppML.threads=4)
if(opt$data == "human_NAc"){
  x <- RcppML::nmf(dat,
                 k=83,
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 1135,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE)
}

if(opt$data == "rat_case_control"){
  x <- RcppML::nmf(dat,
                 k=43,
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 1135,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE)
}

saveRDS(x, file = file.path(res_dir,paste0("nmf_results_", opt$gene_selection_strategy, ".rds")))