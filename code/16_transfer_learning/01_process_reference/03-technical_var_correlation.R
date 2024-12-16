####code for correlating categorical technical variables with NMF patterns
library(here)
library(SpatialExperiment)
library(pheatmap)
library(reshape2)
library(getopt)

spec <- matrix(
    c(
        "gene_selection_strategy", "g", 1, "character", "Choose all genes, or highly deviant genes based on snRNA-seq, or nnSVGs", 
        "data", "d", 1, "character", "Specify input snRNA-seq dataset"
    ),
    byrow = TRUE, ncol = 5
)
opt <- list()
opt$gene_selection_strategy <- "all_genes"
opt$data <- "human_NAc"
opt <- getopt(spec)
print(opt$gene_selection_strategy)

# Read data and create Seurat object
dat_dir <- here::here("processed-data", "12_snRNA")
res_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "RCppML", opt$data)
plot_dir <- here::here("plots", "16_transfer_learning", "01_process_reference", "RCppML", opt$data)


if(opt$data == "human_NAc"){
  sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
}else{
  if(opt$data == "rat_case_control"){
    sce <- readRDS(file = file.path(dat_dir, "NAc_Combo_Integrated.RDS"))
  }else{
        stop("Invalid input data set")
  }
}

sex <- rep("M", dim(sce)[2])
sex[sce$Brain_ID %in% c("Br2720", "Br8325", "Br8492", "Br8667")] <- "F"
colData(sce)$Sex <- sex

x <- readRDS(file = file.path(res_dir,paste0("nmf_results_", opt$gene_selection_strategy, ".rds")))

if(!all.equal(colnames(x@h), colnames(sce))){
  stop("Check the column names of the loadings matrix")
}

## Examining patterns of discrete technical variables
if(opt$data == "human_NAc"){
    data <- as.data.frame(sce$Sort)
    rownames(data) <- colnames(sce)
    colnames(data)<-'Sort'
    onehot_sample <-  dcast(data = data, rownames(data) ~ Sort, length)
    rownames(onehot_sample)<-onehot_sample[,1]
    onehot_sample[,1]<-NULL
    onehot_sample <- onehot_sample[match(rownames(t(x@h)) , rownames(onehot_sample)), ]

    ###correlate with nmf patterns
    pdf(file.path(plot_dir, paste0("nmf_sort_correlation_heatmap_", opt$gene_selection_strategy, ".pdf")))
    pheatmap(cor(t(x@h),onehot_sample), fontsize_row = 5)
    dev.off()

    # Repeat with sample 
    data <- as.data.frame(sce$Sample)
    rownames(data) <- colnames(sce)
    colnames(data)<-'Sample'
    onehot_sample <-  dcast(data = data, rownames(data) ~ Sample, length)
    rownames(onehot_sample)<-onehot_sample[,1]
    onehot_sample[,1]<-NULL
    onehot_sample <- onehot_sample[match(rownames(t(x@h)) , rownames(onehot_sample)), ]


    ###correlate with nmf patterns
    pdf(file.path(plot_dir, paste0("nmf_sample_correlation_heatmap_", opt$gene_selection_strategy, ".pdf")))
    pheatmap(cor(t(x@h),onehot_sample), fontsize_row = 5)
    dev.off()

    # Repeat with Brain ID
    data <- as.data.frame(sce$Brain_ID)
    rownames(data) <- colnames(sce)
    colnames(data)<-'Brain_ID'
    onehot_sample <-  dcast(data = data, rownames(data) ~ Brain_ID, length)
    rownames(onehot_sample)<-onehot_sample[,1]
    onehot_sample[,1]<-NULL
    onehot_sample <- onehot_sample[match(rownames(t(x@h)) , rownames(onehot_sample)), ]

    ###correlate with nmf patterns
    pdf(file.path(plot_dir, paste0("nmf_Brain_ID_correlation_heatmap_", opt$gene_selection_strategy, ".pdf")))
    pheatmap(cor(t(x@h),onehot_sample), fontsize_row = 5)
    dev.off()

    # Repeat with Chromium cDNA date
    data <- as.data.frame(sce$Chromium_cDNA_date)
    rownames(data) <- colnames(sce)
    colnames(data)<-'Chromium_cDNA_date'
    onehot_sample <-  dcast(data = data, rownames(data) ~ Chromium_cDNA_date, length)
    rownames(onehot_sample)<-onehot_sample[,1]
    onehot_sample[,1]<-NULL
    onehot_sample <- onehot_sample[match(rownames(t(x@h)) , rownames(onehot_sample)), ]


    ###correlate with nmf patterns
    pdf(file.path(plot_dir, paste0("nmf_Chromium_cDNA_date_correlation_heatmap_", opt$gene_selection_strategy, ".pdf")))
    pheatmap(cor(t(x@h),onehot_sample), fontsize_row = 5)
    dev.off()

    # Repeat with chromium kit 
    data <- as.data.frame(sce$Chromium_kit)
    rownames(data) <- colnames(sce)
    colnames(data)<-'Chromium_kit'
    onehot_sample <-  dcast(data = data, rownames(data) ~ Chromium_kit, length)
    rownames(onehot_sample)<-onehot_sample[,1]
    onehot_sample[,1]<-NULL
    onehot_sample <- onehot_sample[match(rownames(t(x@h)) , rownames(onehot_sample)), ]

    ###correlate with nmf patterns
    pdf(file.path(plot_dir, paste0("nmf_Chromium_kit_correlation_heatmap_", opt$gene_selection_strategy, ".pdf")))
    pheatmap(cor(t(x@h),onehot_sample), fontsize_row = 5)
    dev.off()

    # Repeat with Sex
    data <- as.data.frame(sce$Sex)
    rownames(data) <- colnames(sce)
    colnames(data)<-'Sex'
    onehot_sample <-  dcast(data = data, rownames(data) ~ Sex, length)
    rownames(onehot_sample)<-onehot_sample[,1]
    onehot_sample[,1]<-NULL
    onehot_sample <- onehot_sample[match(rownames(t(x@h)) , rownames(onehot_sample)), ]

    ###correlate with nmf patterns
    pdf(file.path(plot_dir, paste0("nmf_Sex_correlation_heatmap_", opt$gene_selection_strategy, ".pdf")))
    pheatmap(cor(t(x@h),onehot_sample), fontsize_row = 5)
    dev.off()

    ###for continuous tech vars do this:

    vars<-colData(sce)[,c('detected','sum','subsets_Mito_percent')]

    # fix error cor(t(x@h), vars) : 'y' must be numeric
    vars$detected<-as.numeric(vars$detected)
    vars$sum<-as.numeric(vars$sum)
    vars$mito_percent<-as.numeric(vars$subsets_Mito_percent)
    vars$subsets_Mito_percent<-NULL
    vars <- as.matrix(vars)
    vars <- vars[match(rownames(t(x@h)), rownames(vars)), ]
    pdf(file.path(plot_dir,paste0("nmf_qc_correlation_heatmap_", opt$gene_selection_strategy, ".pdf")))
    pheatmap(cor(t(x@h),vars), fontsize_row = 5)
    dev.off()
}

if(opt$data == "rat_case_control"){
     sce$Stim <- factor(sce$Stim)
     data <- as.data.frame(sce$Stim)
     colnames(data)<-'Stim'
     onehot_sample <-  dcast(data = data, rownames(data) ~ Stim, length)
     rownames(onehot_sample)<-onehot_sample[,1]
     onehot_sample[ ,1] <- NULL
     onehot_sample <- onehot_sample[match(rownames(t(x@h)), rownames(onehot_sample)), ]
    
    pdf(file.path(plot_dir, "nmf_Stim_correlation_heatmap.pdf"))
    pheatmap(cor(t(x@h),onehot_sample), fontsize_row = 5)
    dev.off()

    sce$Sex <- factor(sce$Sex)
    data <- as.data.frame(sce$Sex)
    colnames(data)<-'Sex'
    onehot_sample <-  dcast(data = data, rownames(data) ~ Sex, length)
    rownames(onehot_sample)<-onehot_sample[,1]
    onehot_sample[ ,1] <- NULL
    onehot_sample <- onehot_sample[match(rownames(t(x@h)), rownames(onehot_sample)), ]
    
    pdf(file.path(plot_dir, "nmf_Sex_correlation_heatmap.pdf"))
    pheatmap(cor(t(x@h),onehot_sample), fontsize_row = 5)
    dev.off()

    sce$GEM <- factor(sce$GEM)
    data <- as.data.frame(sce$GEM)
    colnames(data)<-'GEM'
    onehot_sample <-  dcast(data = data, rownames(data) ~ GEM, length)
    rownames(onehot_sample)<-onehot_sample[,1]
    onehot_sample[ ,1] <- NULL
    onehot_sample <- onehot_sample[match(rownames(t(x@h)), rownames(onehot_sample)), ]
    
    pdf(file.path(plot_dir, "nmf_GEM_correlation_heatmap.pdf"))
    pheatmap(cor(t(x@h),onehot_sample), fontsize_row = 5)
    dev.off()

    sce$Dataset <- factor(sce$Dataset)
    data <- as.data.frame(sce$Dataset)
    colnames(data)<-'Dataset'
    onehot_sample <-  dcast(data = data, rownames(data) ~ Dataset, length)
    rownames(onehot_sample)<-onehot_sample[,1]
    onehot_sample[ ,1] <- NULL
    onehot_sample <- onehot_sample[match(rownames(t(x@h)), rownames(onehot_sample)), ]
    
    pdf(file.path(plot_dir, "nmf_Dataset_correlation_heatmap.pdf"))
    pheatmap(cor(t(x@h),onehot_sample), fontsize_row = 5)
    dev.off()

    ###for continuous tech vars do this:

    vars<-sce@meta.data[ ,c("nFeature_RNA", "nCount_RNA", "percent_mito")]

    # fix error cor(t(x@h), vars) : 'y' must be numeric
    vars$nFeature_RNA<-as.numeric(vars$nFeature_RNA)
    vars$nCount_RNA<-as.numeric(vars$nCount_RNA)
    vars$percent_mito<-as.numeric(vars$percent_mito)
    vars <- as.matrix(vars)
    vars <- vars[match(rownames(t(x@h)) ,rownames(vars)), ]

    pdf(file.path(plot_dir,"nmf_qc_correlation_heatmap.pdf"))
    pheatmap(cor(t(x@h),vars), fontsize_row = 5)
    dev.off()
}
