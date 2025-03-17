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
    c( "data", "d", 1, "character", "Specify the input dataset"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

opt <- list()
opt$data <- "human_NAc"
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
  if(opt$data == "rat_case_control_acute"){
    sce <- readRDS(file = file.path(dat_dir, "NAc_Combo_Acute.RDS"))
  }else{
    if(opt$data == "rat_case_control_repeated"){
      sce <- readRDS(file = file.path(dat_dir, "NAc_Combo_Repeated.RDS"))
    }else{
      stop("Invalid input data set")
    }
    
  }
}

x <- readRDS(file = file.path(res_dir,paste0("nmf_results.rds")))

# function for getting top n genes for each pattern
top_genes <- function(W, n=10){
    top_genes <- apply(W, 2, function(x) names(sort(x, decreasing=TRUE)[1:n]))
    return(top_genes)
}
if(opt$data == "human_NAc"){
  geneData <- rowData(sce)
  geneData <- geneData[geneData$gene_id %in% rownames(x@w), ]
  geneData <- geneData[match(rownames(x@w), geneData$gene_id), ]
  rownames(x@w) <- geneData$gene_name
}

# get top 50 genes
top50 <- top_genes(x@w, 50)
write.csv(top50, file = file.path(res_dir, "top50_genes.csv"))

#rat_human_ortholog <- read.delim("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/HOM_AllOrganism.rpt")
#hom_rat <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "rat", ]
#hom_hs <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "human", ]

#top50_human <- data.frame(matrix(NA, nrow = dim(top50)[1], ncol = dim(top50)[2]))
#colnames(top50_human) <- colnames(top50)

#for(i in 1:dim(top50)[2]){
#  print(i)
  #Only run the code if a homolog is found. 
#  if(length(hom_rat[hom_rat$DB.Class.Key %in% hom_hs[hom_hs$Symbol %in% DEGs_df[i,"Gene_symbol"],"DB.Class.Key"],"Symbol"]) == 1){
#    DEGs_df[i,"Rat_Gene_Symbol"] <- hom_rat[hom_rat$DB.Class.Key %in% hom_hs[hom_hs$Symbol %in% DEGs_df[i,"Gene_symbol"],"DB.Class.Key"],"Symbol"]
#  }else{
#    next
#  }
#}

rownames(sce) <- rowData(sce)$gene_name
sce <- sce[rownames(sce) %in% rownames(x@w), ]
sce <- sce[match(rownames(x@w), rownames(sce)), ]
expr <- logcounts(sce)

if(!all.equal(rownames(expr), rownames(x@w))){
  stop("Rownames of expression and NMF loadings do not match")
}

corr_mat <- matrix(NA, nrow = dim(expr)[1], ncol = dim(x@h)[1])
for(i in c(1:dim(x@h)[1])){
  cat(i, "\n")
  corr_mat[ ,i] <- apply(expr, 1, function(igene) cor(igene, x@h[i, ]))
}
rownames(corr_mat) <- rownames(expr)
colnames(corr_mat) <- rownames(x@h)

saveRDS(corr_mat, file.path(res_dir, "gene_corr.rds"))
corr_mat <- reshape2::melt(corr_mat)
corr_mat <- corr_mat[!is.na(corr_mat$value), ]
corr_mat <- corr_mat[abs(corr_mat$value) > 0.2, ]
colnames(corr_mat) <- c("gene", "factor", "correlation")

write.csv(corr_mat, file.path(res_dir, paste0("gene_corr_thresholded_", opt$gene_selection_strategy, ".csv")), row.names= FALSE, quote = FALSE)