library(ggplot2)
library(here)
library(dplyr)
# read data 
dat_dir <- here("processed-data", "19_sLDSC", "10x_visium")
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

# FDR
dat$p_zcore <- pnorm(abs(dat$Coefficient_z.score),lower.tail=F)*2
dat$FDR <- p.adjust(dat$p_zcore,method="fdr")

sig.ldsc <- filter(dat, FDR<.1) 

sig.ldsc$trait.ordered = factor(sig.ldsc$trait, levels=c("Education Years","Intelligence","Neuroticism",  "Bipolar","Schizophrenia","Depression",
                                                         "Smoking initiation","Smoking cessation", "Cigarettes per day", "Drinks per week","Alzheimer Disease", "BMI", "Anorexia"))


palette1= RColorBrewer::brewer.pal(n=7,"RdYlBu")[7:1]
plotDir <- here("plots", "19_sLDSC", "10x_visium")

sig.ldsc$cell.ordered = gsub("_", " ", sig.ldsc$cell)
sig.ldsc$cell.ordered[sig.ldsc$cell.ordered  == "Endothelial Ependymal"] <- "Endothelial/Ependymal"
sig.ldsc$cell.ordered = factor(sig.ldsc$cell.ordered, levels=c("MSN 1","MSN 2", "MSN 3", "D1 islands","Inhibitory", "Excitatory", "Endothelial/Ependymal", "WM"))

pdf(file.path(plotDir, "significant_traits_new.pdf"), width = 7, height = 5)
ggplot(sig.ldsc, aes(x=cell.ordered, y=trait.ordered, size=-log10(FDR), color=Coefficient_z.score))+
  geom_count()+theme_bw()+
  scale_size(range=c(1,4))+
  scale_color_gradientn(colors=c(palette1[1:2],"lightgrey",palette1[6:7]),
                        values=c(0,.2,.4,.8,1))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) + xlab("Spatial Domain") + ylab("Trait") 
dev.off()


