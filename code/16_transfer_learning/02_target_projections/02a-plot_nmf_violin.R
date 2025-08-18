library(SingleCellExperiment)
library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(scater)
library(reshape2)
library(here)

set.seed(123)

opt <- list()
opt$data <- "rat_case_control_acute"
dat_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "preliminary_analysis", opt$data)
sce <- readRDS(file = file.path(dat_dir, "snRNA_seq_NAc.rds"))

# Load the snRNA-seq data
nmf_dir <- here::here('processed-data', '16_transfer_learning')
nmf <- readRDS(file.path(nmf_dir, "01_process_reference", "RCppML", opt$data, paste0("nmf_results.rds")))
loadings <- nmf@h
factors <- nmf@w

loadings <- loadings[ ,match(colnames(sce), colnames(loadings))]
loadings <- t(loadings)

sce@meta.data <- cbind(sce@meta.data, loadings)

if(opt$data == "rat_case_control_acute" | opt$data == "rat_case_control_repeated"){
  data <- as.data.frame(sce$Stim)
  rownames(data) <- colnames(sce)
  colnames(data)<-'Stim'
  onehot_sample <-  dcast(data = data, rownames(data) ~ Stim, length)
  rownames(onehot_sample)<-onehot_sample[,1]
  onehot_sample[,1]<-NULL
  onehot_sample <- onehot_sample[match(rownames(loadings) , rownames(onehot_sample)), ]
  corr_df <- cor(loadings, onehot_sample)
}
if(opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
  data <- as.data.frame(sce$treatment)
  rownames(data) <- colnames(sce)
  colnames(data)<-'Treatment'
  onehot_sample <-  dcast(data = data, rownames(data) ~ Treatment, length)
  rownames(onehot_sample)<-onehot_sample[,1]
  onehot_sample[,1]<-NULL
  onehot_sample <- onehot_sample[match(rownames(loadings) , rownames(onehot_sample)), ]
  corr_df <- cor(loadings, onehot_sample)
}


plotDir <- here::here('plots', '16_transfer_learning', '02_target_projections', opt$data)
sce$CellType_Stim <- paste0(sce$Combo_CellType, "_", sce$Stim)
df <- sce@meta.data

# Extract NMF columns
nmf_cols <- grep("^nmf", colnames(df), value = TRUE)

# Reshape the data to long format
nmf_long <- df %>%
  select(CellType_Stim, all_of(nmf_cols)) %>%
  pivot_longer(cols = all_of(nmf_cols),
               names_to = "NMF_Component",
               values_to = "Weight")

# Optional: remove NAs
nmf_long <- nmf_long %>%
  filter(!is.na(Weight))

pdf(file.path(plotDir, "nmf_cell_type_vlnplots.pdf"), width = 8, height = 3)
for (comp in unique(nmf_long$NMF_Component)) {
  p <- nmf_long %>%
    filter(NMF_Component == comp) %>%
    ggplot(aes(x = CellType_Stim, y = Weight, fill = CellType_Stim)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.5) +
    theme_minimal() +
    labs(title = paste("Violin Plot for", comp),
         x = "CellType_Stim",
         y = "NMF Weight") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  print(p)
}
dev.off()