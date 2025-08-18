library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(here)
set.seed(123)

# read data 
dat_dir <- here("processed-data", "19_sLDSC", "human_NAc_NMF")
dat <- read.table(file.path(dat_dir, "ldsc_results.txt"),as.is=T,header=T,sep="\t")

traits <- c(
"ADHD",
"Alzheimer Disease3",
 "Anorexia",
 "Autism",
 "BMI",
 "Bipolar Disorder2",
 "GSCAN_AgeSmk",
 "GSCAN_CigDay",
 "GSCAN_DrnkWk",
 "GSCAN_SmkCes",
 "GSCAN_SmkInit",
 "Height",
 "Type_2_Diabetes",
 "mdd2019edinburgh",
 "epilepsy",
 "Intelligence",
 "Parkinson Disease",
 "Schizophrenia_PGC3",
 "Education Years",
 "Neuroticism"
)
idx <- is.element(dat$trait,traits)
dat <- dat[idx,]

# rename traits
dat$trait[dat$trait=="GSCAN_AgeSmk"] <- "Age of smoking"
dat$trait[dat$trait=="GSCAN_CigDay"] <- "Cigarettes per day"
dat$trait[dat$trait=="GSCAN_DrnkWk"] <- "Drinks per week"
dat$trait[dat$trait=="GSCAN_SmkCes"] <- "Smoking cessation"
dat$trait[dat$trait=="GSCAN_SmkInit"] <- "Smoking initiation"
dat$trait[dat$trait=="Schizophrenia_PGC3"] <- "Schizophrenia"
dat$trait[dat$trait=="Bipolar Disorder2"] <- "Bipolar"
dat$trait[dat$trait=="epilepsy"] <- "Epilepsy"
dat$trait[dat$trait=="Type_2_Diabetes"] <- "Type 2 Diabetes"
dat$trait[dat$trait=="Alzheimer Disease3"] <- "Alzheimer Disease"
dat$trait[dat$trait=="mdd2019edinburgh"] <- "Depression"

dat$p_zcore <- pnorm(abs(dat$Coefficient_z.score),lower.tail=F)*2
dat$FDR <- p.adjust(dat$p_zcore,method="fdr")

# Read in the NMF inputs
res_dir <- here("processed-data", "16_transfer_learning", "02_target_projections", "human_NAc")
spe <- readRDS(file = file.path(res_dir, paste0("spe_NMF.rds")))

non0.spots = colSums(as.matrix(colData(spe)[,paste0("nmf",1:66)])>0)
remove.nmf = names(non0.spots[non0.spots<200])
remove.nmf = c(remove.nmf, "nmf26", "nmf45")

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

sig.ldsc = filter(dat, cell %in% nmf.ordered.keep, FDR<.05) 
sig.ldsc$trait.ordered = factor(sig.ldsc$trait,
 levels=c("BMI","Height","Education Years","Drinks per week",
          "Smoking initiation","Smoking cessation",
           "Intelligence","Neuroticism",
           "Epilepsy","Schizophrenia",
           "Bipolar","Depression",
           "Alzheimer Disease"))
sig.ldsc$nmf_f= factor(sig.ldsc$cell, levels=nmf.ordered.keep)

palette1= RColorBrewer::brewer.pal(n=7,"RdYlBu")[7:1]
plotDir <- here("plots", "19_sLDSC", "human_NAc_NMF")

pdf(file.path(plotDir, "ldsc_results.pdf"), width = 10, height = 4)
ggplot(sig.ldsc, aes(x=nmf_f, y=trait.ordered, size=-log10(FDR), color=Coefficient_z.score))+
  geom_count()+theme_bw()+
  scale_size(range=c(1,4))+
  #scale_color_gradient2(low=palette1[1], mid="lightgrey", high=palette1[7])+
  scale_color_gradientn(colors=c(palette1[1:2],"lightgrey",palette1[6:7]),
                        values=c(0,.2,.4,.8,1))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) + xlab("") + ylab("")
dev.off()                                       

MSN_select <- c("nmf38", "nmf10", "nmf3", "nmf7", "nmf39", "nmf4", "nmf25", "nmf8", 
"nmf9", "nmf15", "nmf31", "nmf50", "nmf65", "nmf56")
sig.ldsc_MSN <- sig.ldsc[sig.ldsc$nmf %in% MSN_select, ]
sig.ldsc_MSN$nmf_f <- factor(sig.ldsc_MSN$nmf, levels= MSN_select)
pdf(file.path(plotDir, "ldsc_results_MSN.pdf"), width = 5, height = 4)
ggplot(sig.ldsc_MSN, aes(x=nmf_f, y=trait.ordered, size=-log10(FDR), color=Coefficient_z.score))+
  geom_count()+theme_bw()+
  scale_size(range=c(1,4))+
  #scale_color_gradient2(low=palette1[1], mid="lightgrey", high=palette1[7])+
  scale_color_gradientn(colors=c(palette1[1:2],"lightgrey",palette1[6:7]),
                        values=c(0,.2,.4,.8,1))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) + xlab("") + ylab("")
dev.off()   

D1_select <- c("nmf34", "nmf35", "nmf44", "nmf8", "nmf9", "nmf31", "nmf64")
sig.ldsc_D1 <- sig.ldsc[sig.ldsc$nmf %in% D1_select, ]
sig.ldsc_D1$nmf_f <- factor(sig.ldsc_D1$nmf, levels= D1_select)
pdf(file.path(plotDir, "ldsc_results_D1.pdf"), width = 5, height = 4)
ggplot(sig.ldsc_D1, aes(x=nmf_f, y=trait.ordered, size=-log10(FDR), color=Coefficient_z.score))+
  geom_count()+theme_bw()+
  scale_size(range=c(1,4))+
  #scale_color_gradient2(low=palette1[1], mid="lightgrey", high=palette1[7])+
  scale_color_gradientn(colors=c(palette1[1:2],"lightgrey",palette1[6:7]),
                        values=c(0,.2,.4,.8,1))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) + xlab("") + ylab("")
dev.off()   