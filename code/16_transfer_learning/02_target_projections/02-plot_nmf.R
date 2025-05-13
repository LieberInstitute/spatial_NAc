library(SingleCellExperiment)
library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(scater)
library(reshape2)
library(here)

set.seed(123)

opt <- list()
opt$data <- "rat_case_control_repeated"

dat_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "preliminary_analysis", opt$data)
sce <- readRDS(file = file.path(dat_dir, "snRNA_seq_NAc.rds"))

# Load the snRNA-seq data
nmf_dir <- here::here('processed-data', '16_transfer_learning')
nmf <- readRDS(file.path(nmf_dir, "01_process_reference", "RCppML", opt$data, paste0("nmf_results.rds")))
loadings <- nmf@h
factors <- nmf@w

if(opt$data == "human_NAc"){
  sex <- rep("M", dim(sce)[2])
  sex[sce$Brain_ID %in% c("Br2720", "Br8325", "Br8492", "Br8667")] <- "F"
  sce$Sex <- sex
}

loadings <- loadings[ ,match(colnames(sce), colnames(loadings))]
loadings <- t(loadings)

sce@meta.data <- cbind(sce@meta.data, loadings)

if(opt$data == "human_NAc"){
  data <- as.data.frame(sce$Sex)
  rownames(data) <- colnames(sce)
  colnames(data)<-'Sex'
  onehot_sample <-  dcast(data = data, rownames(data) ~ Sex, length)
  rownames(onehot_sample)<-onehot_sample[,1]
  onehot_sample[,1]<-NULL
  onehot_sample <- onehot_sample[match(rownames(loadings) , rownames(onehot_sample)), ]
  corr_df <- cor(loadings, onehot_sample)
}
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



if(opt$data == "human_NAc"){
  corr_df <- corr_df[abs(corr_df[ ,1]) > 0.3, ]
}
if(opt$data == "rat_case_control_acute" | opt$data == "rat_case_control_repeated"){
  corr_df <- corr_df[abs(corr_df[ ,1]) > 0.1, ]
}
if(opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
  corr_df <- corr_df[abs(corr_df[ ,1]) > 0.3, ]
}

plotDir <- here::here('plots', '16_transfer_learning', '02_target_projections', opt$data)
#baseline non-zero
if(opt$data == "human_NAc"){
  nFactors <- 66
}
if(opt$data == "rat_case_control_acute"){
  nFactors <- 43
}
if(opt$data == "rat_case_control_repeated"){
  nFactors <- 40
}
if(opt$data == "rat_case_control_morphine_acute"){
  nFactors <- 46
}
if(opt$data == "rat_case_control_morphine_repeated"){
  nFactors <- 36
}

meta_df <- sce@meta.data
non0.nuc = colSums(as.matrix(meta_df[,paste0("nmf",1:nFactors)])>0)

pdf(file.path(plotDir, "weighted_nuclei_per_NMF.pdf"), width = 5, height = 4)
plot(ecdf(log10(non0.nuc)), xlim=c(0,5), xlab="# non-zero weighted nuclei per NMF", main = "CDF of log10(#nuclei)")
dev.off()
# Load the spe object 
spe <- readRDS(file.path(nmf_dir, "02_target_projections", opt$data, paste0("spe_NMF.rds")))
clusters_resFile <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters", "precast_clusters.csv")
spe[["domain"]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(clusters_resFile), by = 'key') |>
    pull(cluster) |>
    as.factor()

non0.spots = colSums(as.matrix(colData(spe)[,paste0("nmf",1:nFactors)])>0)
plot.spots = ifelse(na.exclude(non0.spots < 200), "red", "black")
pdf(file.path(plotDir, "weighted_spots_per_NMF.pdf"), width = 5, height = 4)
plot(ecdf(log10(non0.spots)), xlim=c(0,5), xlab="# non-zero weighted spots per NMF",
     col=sort(plot.spots, decreasing=T), 
     main="NMF with < 200 spots removed")
dev.off()

remove.nmf = names(non0.spots[non0.spots<200])

##############
## dotplot
if(opt$data == "rat_case_control_acute" | opt$data == "rat_case_control_repeated"){
  sce$cellType_Stim <- paste(sce$Combo_CellType, sce$Stim, sep = " ")
}

if(opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
  sce$cellType_Stim <- paste(sce$CellType, sce$treatment, sep = " ")
}
meta_df <- sce@meta.data

if(opt$data == "human_NAc"){
  cellType_col <- "CellType"
}
if(opt$data == "rat_case_control_acute" | opt$data == "rat_case_control_repeated"){
  cellType_col <- "cellType_Stim"
}
if(opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
  cellType_col <- "cellType_Stim"
}
seed1 = as.matrix(meta_df[,paste0("nmf",1:nFactors)])
seed1 = seed1>0
d1 = cbind.data.frame(cell.class=as.data.frame(meta_df)[,cellType_col],
                      seed1) %>% 
  group_by(cell.class) %>% add_tally(name="total") %>%
    group_by(cell.class, total) %>%
  summarise_at(paste0("nmf",1:nFactors), sum) %>%
  tidyr::pivot_longer(paste0("nmf",1:nFactors), values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)

seed2 = as.matrix(meta_df[,paste0("nmf",1:nFactors)])
seed2 = apply(seed2, 2, scale)
d2 = cbind.data.frame(cell.class=as.data.frame(meta_df)[,cellType_col],
                      seed2) %>% 
  group_by(cell.class) %>%
  summarise_at(paste0("nmf",1:nFactors), mean) %>% 
  tidyr::pivot_longer(paste0("nmf",1:nFactors), values_to="scaled.avg", names_to="nmf")

if(opt$data == "human_NAc"){
  spf.ordered <- c("DRD1_MSN_C", "DRD1_MSN_A", "DRD2_MSN_A", "DRD2_MSN_B", "DRD1_MSN_B", "DRD1_MSN_D", 
  "Inh_A", "Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F",
  "Excitatory", "Astrocyte_A", "Astrocyte_B", "Oligo", "OPC", 
  "Microglia", "Endothelial", "Ependymal", "Neuron_Ambig")
}

if(opt$data == "rat_case_control_acute" | opt$data == "rat_case_control_repeated"){
  spf.ordered <- c("Drd1-MSN-1 Cocaine","Drd1-MSN-1 Saline",
   "Drd2-MSN-1 Cocaine",  "Drd2-MSN-1 Saline","Drd1-MSN-2 Cocaine",
   "Drd1-MSN-2 Saline","Drd2-MSN-2 Cocaine", "Drd2-MSN-2 Saline",
  "Drd3-MSN Cocaine", "Drd3-MSN Saline", "Grm8-MSN Cocaine", "Grm8-MSN Saline",
  "Glutamatergic Cocaine","Glutamatergic Saline", "GABAergic Cocaine","GABAergic Saline",
   "Chat-Interneuron Cocaine", "Chat-Interneuron Saline","Sst-Interneuron Cocaine","Sst-Interneuron Saline",
    "Pvalb-Interneuron Cocaine","Pvalb-Interneuron Saline",
    "Astrocyte Cocaine","Astrocyte Saline", 
    "Polydendrocyte Cocaine","Polydendrocyte Saline", "Olig-1 Cocaine","Olig-1 Saline", 
     "Microglia Cocaine",  "Microglia Saline", "Mural Cocaine", "Mural Saline"
   )
}

if(opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
  spf.ordered <- c("MSN1 morphine", "MSN1 saline", "MSN2 morphine", "MSN2 saline", 
   "MSN3 morphine", "MSN3 saline",
   "Drd1 MSN1 morphine", "Drd1 MSN1 saline",
   "Drd2 MSN1 morphine","Drd2 MSN1 saline",
   "Drd1 MSN2 morphine","Drd1 MSN2 saline",
   "Drd2 MSN2 morphine","Drd2 MSN2 saline",
   "Drd1 MSN3 morphine", "Drd1 MSN3 saline",
   "Drd2 MSN3 morphine", "Drd2 MSN3 saline",
   "Drd2 MSN4 morphine", "Drd2 MSN4 saline", 
   "Grm8 MSN morphine", "Grm8 MSN saline",
   "Cholinergic morphine", "Cholinergic saline",
   "Sst Inter morphine","Sst Inter saline", 
   "Pvalb Inter morphine","Pvalb Inter saline",
   "Interneuron1 morphine", "Interneuron1 saline", 
   "Interneuron2 morphine", "Interneuron2 saline",
   "Interneuron3 morphine", "Interneuron3 saline",
   "Astro1 morphine", "Astro1 saline",
   "Astro2 morphine","Astro2 saline",
   "Astro3 morphine","Astro3 saline",
  "Oligo1 morphine", "Oligo1 saline", "Oligo2 morphine","Oligo2 saline",
   "Oligo3 morphine","Oligo3 saline",
   "Oligo4 morphine","Oligo4 saline",
   "OPC morphine","OPC saline",
   "Micro morphine","Micro saline",
   "Endo morphine","Endo saline"
   )
}

dot.df = left_join(d1[,c("cell.class","nmf","prop")], 
                   d2[,c("cell.class","nmf","scaled.avg")]) %>%
  mutate(cell.class=factor(cell.class, 
                                     levels=spf.ordered))


if(opt$data == "human_NAc"){
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
}

if(opt$data == "rat_case_control_acute"){
  nmf.ordered = c( "nmf1", "nmf5", "nmf9", "nmf36","nmf39", "nmf32",
                   "nmf2", "nmf3", "nmf23", 
                   "nmf4", "nmf6", # General neuron
                   "nmf8", "nmf16" , # Drd1-MSN-1
                   "nmf12" , # Drd2-MSN-1
                   "nmf19", "nmf20" , # Drd1-MSN-2
                   "nmf30" , # Drd2-MSN-2
                   "nmf29", # Drd3 MSN
                   "nmf24", "nmf28", "nmf31", #Grm8-MSN
                   "nmf37" ,# Glutamergic
                   "nmf17",# GABAergic
                   "nmf18", # Chat Interneuron
                   "nmf35" , # SST interneuron
                   "nmf33", "nmf38" , #Pvalb Interneuron
                   "nmf13", "nmf27","nmf42" , # Astrocyte
                   "nmf7","nmf41","nmf40", # Polydendrocyte
                   "nmf10", "nmf11", "nmf14", "nmf15", "nmf21", "nmf22", "nmf26" , # Olig-1
                   "nmf25", "nmf34" , # Microglia
                   "nmf43" #Mural
                   )
}

if(opt$data == "rat_case_control_repeated"){
   nmf.ordered = c("nmf1", "nmf2", "nmf4","nmf22","nmf30","nmf33","nmf35","nmf37", "nmf38","nmf39", # General
                   "nmf11", "nmf3", "nmf5", "nmf6", # General neuron
                   "nmf17", "nmf7", # General MSN
                   "nmf19", #Drd1-MSN-1
                   "nmf14" , # Drd2-MSN-1
                   "nmf15", # Drd1-MSN-2
                   "nmf28" , # Drd2-MSN-2
                   "nmf34", # Drd3 MSN
                   "nmf21", "nmf24","nmf31", #Grm8-MSN
                   "nmf29" ,# Glutamergic
                   "nmf13",# GABAergic
                   "nmf9", # Chat Interneuron
                   "nmf36" , # SST interneuron
                   "nmf32" , #Pvalb Interneuron
                   "nmf10", "nmf26" , # Astrocyte
                   "nmf16","nmf40", # Polydendrocyte
                   "nmf8","nmf12","nmf18","nmf20","nmf23", # Olig-1
                   "nmf25" , # Microglia
                   "nmf27" #Mural
                   )
}
if(opt$data == "rat_case_control_morphine_acute"){
  nmf.ordered = c("nmf2", "nmf5", "nmf6", "nmf7","nmf13","nmf14","nmf18","nmf22","nmf31","nmf36","nmf39","nmf43", "nmf44",  # General
                "nmf8","nmf10","nmf29",  # General neuron
                "nmf1","nmf4","nmf25",  # General MSN
                "nmf24", #MSN 2
                "nmf23", #MSN 3
                "nmf21","nmf12","nmf16", # Drd1 MSN 1
                "nmf30", # Drd2 MSN 1
                "nmf15", #Drd1 MSN 2
                "nmf11","nmf17", # Drd2 MSN 2
                "nmf20", #Drd1 MSN 3
                "nmf34", # Grm8 MSN
                "nmf40", # Cholinergic
                "nmf37", # sst inter
                "nmf32" , #Pvalb Interneuron
                "nmf38", #Interneuron 1
                "nmf33", #Interneuron 2
                "nmf35","nmf3", "nmf27", #Astro
                "nmf41", "nmf42", "nmf9", "nmf26", "nmf19", #Oligo
                "nmf28", #OPC
                "nmf45", #Micro
                "nmf46" #Endo
  )
}
if(opt$data == "rat_case_control_morphine_repeated"){
  nmf.ordered = c("nmf1","nmf10","nmf12", "nmf13","nmf2","nmf26", # General
                  "nmf21","nmf8", # General neuron  
                  "nmf3","nmf5", # General MSN
                  "nmf22", # MSN2
                  "nmf20", #MSN3
                  "nmf14", # Drd1 MSN 1 
                  "nmf6", #Drd1 MSN 2
                  "nmf24", # Drd2 MSN 1
                  "nmf7", # Drd2 MSN 2
                  "nmf17", #Drd2 MSN 2 
                  "nmf16", # Drd1 MSN 3
                  "nmf11", # Drd2 MSN 3
                  "nmf28", #Grm8 MSN
                  "nmf32", #Cholinergic
                  "nmf31", #Sst Inter
                  "nmf27", #Pvalb interneuron
                  "nmf19", # Interneuron 1
                  "nmf25", #Interneuron
                  "nmf29", #Astro 1
                  "nmf4", #Astro 2
                  "nmf30", #Astro 3
                  "nmf33","nmf34","nmf9", #Oligo 1
                  "nmf18", # Oligo 2
                  "nmf15", #Oligo 4
                  "nmf23", #OPC
                  "nmf35", #Micro
                  "nmf36" # Endo
  )
}
             
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


#spe dotplot
spe <- spe[ ,!is.na(spe[["domain"]])]
seed3 = as.matrix(colData(spe)[,c(rownames(corr_df),nmf.ordered.remove, nmf.ordered.keep)])
seed3 = seed3>0
d3 = cbind.data.frame(domain=as.data.frame(colData(spe))[,"domain"],
                      seed3) %>% 
  group_by(domain) %>% add_tally(name="total") %>%
  group_by(domain, total) %>%
  summarise_at(c(rownames(corr_df),nmf.ordered.remove, nmf.ordered.keep), sum) %>%
  tidyr::pivot_longer(c(rownames(corr_df),nmf.ordered.remove, nmf.ordered.keep), values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)

seed4 = as.matrix(colData(spe)[,c(rownames(corr_df),nmf.ordered.remove, nmf.ordered.keep)])
seed4 = apply(seed4, 2, scale)
d4 = cbind.data.frame(domain=as.data.frame(colData(spe))[,"domain"],
                      seed4) %>% 
  group_by(domain) %>%
  summarise_at(c(rownames(corr_df),nmf.ordered.remove, nmf.ordered.keep), mean) %>% 
  tidyr::pivot_longer(c(rownames(corr_df),nmf.ordered.remove, nmf.ordered.keep), values_to="scaled.avg", names_to="nmf")

dot.df2 = left_join(d3[,c("domain","nmf","prop")], 
                    d4[,c("domain","nmf","scaled.avg")])

domain.ordered <- c("MSN 1", "MSN 2", "MSN 3", "D1 islands", "Inhibitory", "Excitatory", 
  "WM", "Endothelial/Ependymal")

dot.df2$domain = factor(dot.df2$domain, levels = domain.ordered)
dot.df2$nmf_f = factor(dot.df2$nmf, levels=c(rownames(corr_df),nmf.ordered.remove, nmf.ordered.keep))

pdf(file.path(plotDir, "NMF_by_spatial_domains.pdf"), width = 11, height = 4)
ggplot(dot.df2, aes(x=nmf_f, y=domain, size=prop, color=scaled.avg))+
  geom_count()+theme_bw()+
  scale_size(range=c(0,3))+scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) + xlab("") + ylab("")
dev.off()

if(opt$data == "human_NAc"){
  dot.df3 <- dot.df2[dot.df2$domain %in% c("MSN 1", "MSN 2", "MSN 3"), ]
dot.df3 <- dot.df3[dot.df3$prop > 0.2 & dot.df3$scaled.avg > 0.2, ]
unique_nmf <- unique(dot.df3$nmf)
select_nmf <- c("nmf38", "nmf10", "nmf3", "nmf7", "nmf39", "nmf4", "nmf25")
unique_nmf_ordered <- c(select_nmf, unique_nmf[!unique_nmf %in% select_nmf])

dot.df4 <- dot.df2[dot.df2$nmf %in% unique_nmf_ordered, ]
dot.df4$nmf_f <- factor(dot.df4$nmf, levels = unique_nmf_ordered)

pdf(file.path(plotDir, "NMF_MSN_domains.pdf"), width = 6, height = 4)
ggplot(dot.df4, aes(x=nmf_f, y=domain, size=prop, color=scaled.avg))+
  geom_count()+theme_bw()+
  scale_size(range=c(0,3))+scale_color_viridis_c(option="F", direction=-1, limits = c(-1.6, 3))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) + xlab("") + ylab("")
dev.off()

dot.df5 <- dot.df2[dot.df2$domain %in% c("D1 islands"), ]
dot.df5 <- dot.df5[dot.df5$prop > 0.2 & dot.df5$scaled.avg > 0.2, ]
unique_nmf <- unique(dot.df5$nmf)
select_nmf <- c("nmf34", "nmf35", "nmf44")
unique_nmf_ordered <- c(select_nmf, unique_nmf[!unique_nmf %in% select_nmf])
dot.df6 <- dot.df2[dot.df2$nmf %in% unique_nmf_ordered, ]
dot.df6$nmf_f <- factor(dot.df6$nmf, levels = unique_nmf_ordered)
pdf(file.path(plotDir, "NMF_D1_islands.pdf"), width = 5, height = 4)
ggplot(dot.df6, aes(x=nmf_f, y=domain, size=prop, color=scaled.avg))+
  geom_count()+theme_bw()+
  scale_size(range=c(0,3))+scale_color_viridis_c(option="F", direction=-1, limits = c(-1.6, 3.1))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) + xlab("") + ylab("")
dev.off()
}
