rm(list = ls())
library(SingleCellExperiment)
library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(scater)
library(reshape2)
library(here)

set.seed(123)

opt <- list()
opt$data <- "rat_case_control_cocaine_repeated"
res_file <- here::here("processed-data", "16_transfer_learning", "02_target_projections", "projection_of_human_factors", paste0(opt$data, ".rds"))
sce <- readRDS(res_file)

if(opt$data == "rat_case_control_cocaine_acute" | opt$data == "rat_case_control_cocaine_repeated"){
  sce$cellType_Stim <- paste(sce$Combo_CellType, sce$Stim, sep = " ")
  df <- table(sce$Combo_CellType, sce$Stim)
  df <- reshape2::melt(df)
}

if(opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
  sce$cellType_Stim <- paste(sce$CellType, sce$treatment, sep = " ")
  df <- table(sce$CellType, sce$treatment)
  df <- reshape2::melt(df)
}

plotDir <- here::here("plots", "16_transfer_learning", "02_target_projections", "projection_human_factors", opt$data)
pdf(file.path(plotDir, "cells_per_condition.pdf"), width = 5, height = 6)
ggplot(df, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity", position = "stack") + xlab("") + ylab("") + theme_classic() + coord_flip()
dev.off()


meta_df <- sce@meta.data
projections <- meta_df[ ,grep("nmf", colnames(meta_df))]
if(opt$data == "rat_case_control_cocaine_acute" | opt$data == "rat_case_control_cocaine_repeated"){
  data <- as.data.frame(sce$Stim)
  rownames(data) <- colnames(sce)
  colnames(data)<-'Stim'
  onehot_sample <-  dcast(data = data, rownames(data) ~ Stim, length)
  rownames(onehot_sample)<-onehot_sample[,1]
  onehot_sample[,1]<-NULL
  onehot_sample <- onehot_sample[match(rownames(projections) , rownames(onehot_sample)), ]
  corr_df <- cor(projections, onehot_sample)
}
if(opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
  data <- as.data.frame(sce$treatment)
  rownames(data) <- colnames(sce)
  colnames(data)<-'Treatment'
  onehot_sample <-  dcast(data = data, rownames(data) ~ Treatment, length)
  rownames(onehot_sample)<-onehot_sample[,1]
  onehot_sample[,1]<-NULL
  onehot_sample <- onehot_sample[match(rownames(projections) , rownames(onehot_sample)), ]
  corr_df <- cor(projections, onehot_sample)
}


if(opt$data == "rat_case_control_cocaine_acute" | opt$data == "rat_case_control_cocaine_repeated"){
  corr_df <- corr_df[abs(corr_df[ ,1]) > 0.1, ,drop = FALSE]
}
if(opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
  corr_df <- corr_df[abs(corr_df[ ,1]) > 0.1, ,drop = FALSE]
}


non0.cells = colSums(as.matrix(meta_df[,paste0("nmf",1:66)])>0)
plot.cells = ifelse(na.exclude(non0.cells < 2000), "red", "black")
pdf(file.path(plotDir, "weighted_cells_per_NMF.pdf"), width = 5, height = 4)
plot(ecdf(log10(non0.cells)), xlim=c(0,5), xlab="# non-zero weighted spots per NMF",
     col=sort(plot.cells, decreasing=T), 
     main="NMF with < 2000 cells removed")
dev.off()

remove.nmf = names(non0.cells[non0.cells<2000])
if(opt$data == "rat_case_control_cocaine_acute" | opt$data == "rat_case_control_cocaine_repeated"){
  cellType_col <- "cellType_Stim"
}
if(opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
  cellType_col <- "cellType_Stim"
}

seed1 = as.matrix(meta_df[,paste0("nmf",1:66)])
seed1 = seed1>0
d1 = cbind.data.frame(cell.class=as.data.frame(meta_df)[,cellType_col],
                      seed1) %>% 
  group_by(cell.class) %>% add_tally(name="total") %>%
    group_by(cell.class, total) %>%
  summarise_at(paste0("nmf",1:66), sum) %>%
  tidyr::pivot_longer(paste0("nmf",1:66), values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)

seed2 = as.matrix(meta_df[,paste0("nmf",1:66)])
seed2 = apply(seed2, 2, scale)
d2 = cbind.data.frame(cell.class=as.data.frame(meta_df)[,cellType_col],
                      seed2) %>% 
  group_by(cell.class) %>%
  summarise_at(paste0("nmf",1:66), mean) %>% 
  tidyr::pivot_longer(paste0("nmf",1:66), values_to="scaled.avg", names_to="nmf")

if(opt$data == "rat_case_control_cocaine_acute" | opt$data == "rat_case_control_cocaine_repeated"){
  spf.ordered <- c("Drd1-MSN-1 Cocaine", "Drd1-MSN-1 Saline",
                   "Drd2-MSN-1 Cocaine", "Drd2-MSN-1 Saline",
                   "Drd1-MSN-2 Cocaine", "Drd1-MSN-2 Saline",
                   "Drd2-MSN-2 Cocaine", "Drd2-MSN-2 Saline",
                   "Grm8-MSN Cocaine", "Grm8-MSN Saline",
                   "Drd3-MSN Cocaine", "Drd3-MSN Saline",
                   "GABAergic Cocaine", "GABAergic Saline",
                   "Chat-Interneuron Cocaine", "Chat-Interneuron Saline",
                   "Pvalb-Interneuron Cocaine", "Pvalb-Interneuron Saline",
                   "Sst-Interneuron Cocaine", "Sst-Interneuron Saline",
                   "Glutamatergic Cocaine", "Glutamatergic Saline", 
                   "Astrocyte Cocaine", "Astrocyte Saline",
                   "Olig-1 Cocaine", "Olig-1 Saline",
                   "Polydendrocyte Cocaine", "Polydendrocyte Saline",
                   "Microglia Cocaine", "Microglia Saline",
                   "Mural Cocaine", "Mural Saline"
   )
}
if(opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
  spf.ordered <- c("MSN1 morphine", "MSN1 saline", "MSN2 morphine", "MSN2 saline", "MSN3 morphine", "MSN3 saline",
   "Drd1 MSN1 morphine","Drd1 MSN1 saline", "Drd2 MSN1 morphine","Drd2 MSN1 saline",
   "Drd1 MSN2 morphine", "Drd1 MSN2 saline", "Drd2 MSN2 morphine","Drd2 MSN2 saline",
    "Drd1 MSN3 morphine", "Drd1 MSN3 saline", "Drd2 MSN3 morphine","Drd2 MSN3 saline",
    "Drd2 MSN4 morphine", "Drd2 MSN4 saline", "Grm8 MSN morphine", "Grm8 MSN saline",
    "Pvalb Inter morphine","Pvalb Inter saline", "Interneuron2 morphine", "Interneuron2 saline",
     "Interneuron3 morphine","Interneuron3 saline","Interneuron1 morphine","Interneuron1 saline",
      "Cholinergic morphine","Cholinergic saline", "Sst Inter morphine","Sst Inter saline",
  "Astro2 morphine","Astro2 saline", "Astro1 morphine", "Astro1 saline",
  "Astro3 morphine","Astro3 saline", "Oligo2 morphine","Oligo2 saline", "Oligo3 morphine",
  "Oligo3 saline", "Oligo4 morphine", "Oligo4 saline", "Oligo1 morphine","Oligo1 saline", 
  "OPC morphine", "OPC saline","Micro morphine","Micro saline", "Endo morphine", "Endo saline"
   )
}
dot.df = left_join(d1[,c("cell.class","nmf","prop")], 
                   d2[,c("cell.class","nmf","scaled.avg")]) %>%
  mutate(cell.class=factor(cell.class, 
                                     levels=spf.ordered))

# Specify the order of the NMF
nmf.ordered = c("nmf45","nmf1","nmf21","nmf48","nmf12", "nmf14","nmf16","nmf23","nmf28","nmf30", "nmf40", "nmf41", "nmf43","nmf51", "nmf52","nmf56","nmf66", # General
                "nmf26","nmf5","nmf8", "nmf9","nmf15","nmf25","nmf31","nmf38",   # General neuron
                "nmf18", # General MSN
                "nmf24", # DRD1 MSN C
                "nmf6", "nmf11", "nmf7", "nmf2", "nmf3", "nmf4","nmf10",  # DRD1/DRD2 MSN A
                "nmf37", # DRD2 MSN B
                "nmf34", "nmf36", "nmf35", # DRD1 MSN B
                "nmf44", # DRD1 MSN D
                "nmf32", # Inhib A
                "nmf57", # Inhib B
                "nmf42", # Inhib C
                "nmf59", # Inhib D
                "nmf58", # Inhib E
                "nmf61", #Inhib F
                "nmf49", # Excitatory
                "nmf20", "nmf33", "nmf39", "nmf50","nmf60", "nmf54", # Astrocyte A/B
                "nmf13", "nmf19", "nmf22", "nmf27", "nmf46", "nmf47", #Oligo
                "nmf17", #OPC
                "nmf29", "nmf55", "nmf62", "nmf64",  # Microglia
                "nmf63", "nmf65", # Endothelial
                "nmf53" # Ependymal
                )

nmf.ordered.keep = setdiff(nmf.ordered, remove.nmf)
nmf.ordered.remove = intersect(nmf.ordered, remove.nmf)
sex_specific_factors <- rownames(corr_df)
nmf.ordered.keep <- nmf.ordered.keep[!nmf.ordered.keep %in% sex_specific_factors]
nmf.ordered.remove <- nmf.ordered.remove[!nmf.ordered.remove %in% sex_specific_factors]

dot.df$nmf_f = factor(dot.df$nmf, levels=c(rownames(corr_df), nmf.ordered.remove, nmf.ordered.keep))
pdf(file.path(plotDir, "NMF_by_cell_types.pdf"), width = 10, height = 8)
ggplot(dot.df, aes(x=nmf_f, y=cell.class, size=prop, color=scaled.avg))+
  geom_count()+theme_bw()+
  scale_size(range=c(0,3))+scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) + xlab("") + ylab("")
dev.off()
