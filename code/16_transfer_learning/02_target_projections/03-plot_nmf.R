library(SingleCellExperiment)
library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(scater)

set.seed(123)

opt <- list()
opt$data <- "human_NAc"
opt$gene_selection_strategy <- "all_genes"

# Load the snRNA-seq data
nmf_dir <- here::here('processed-data', '16_transfer_learning')
sc_dir <- here::here("processed-data", "12_snRNA")
sce <- readRDS(file = file.path(sc_dir, "sce_CellType_noresiduals.Rds"))
nmf <- readRDS(file.path(nmf_dir, "01_process_reference", "RCppML", opt$data, paste0("nmf_results_", opt$gene_selection_strategy, ".rds")))
loadings <- nmf@h
factors <- nmf@w

sex <- rep("M", dim(sce)[2])
sex[sce$Brain_ID %in% c("Br2720", "Br8325", "Br8492", "Br8667")] <- "F"
colData(sce)$Sex <- sex

loadings <- loadings[ ,match(colnames(sce), colnames(loadings))]
loadings <- t(loadings)
colData(sce) <- cbind(colData(sce), loadings)

data <- as.data.frame(sce$Sex)
rownames(data) <- colnames(sce)
colnames(data)<-'Sex'
onehot_sample <-  dcast(data = data, rownames(data) ~ Sex, length)
rownames(onehot_sample)<-onehot_sample[,1]
onehot_sample[,1]<-NULL
onehot_sample <- onehot_sample[match(rownames(loadings) , rownames(onehot_sample)), ]
corr_df <- cor(loadings, onehot_sample)
corr_df <- corr_df[abs(corr_df[ ,1]) > 0.3, ]

plotDir <- here::here('plots', '16_transfer_learning', '02_target_projections', opt$data)
#baseline non-zero
nFactors <- 83
non0.nuc = colSums(as.matrix(colData(sce)[,paste0("nmf",1:nFactors)])>0)
pdf(file.path(plotDir, "weighted_nuclei_per_NMF.pdf"), width = 5, height = 4)
plot(ecdf(log10(non0.nuc)), xlim=c(0,5), xlab="# non-zero weighted nuclei per NMF", main = "CDF of log10(#nuclei)")
dev.off()
# Load the spe object 
spe <- readRDS(file.path(nmf_dir, "02_target_projections", opt$data, paste0("spe_NMF_", opt$gene_selection_strategy, ".rds")))
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
seed1 = as.matrix(colData(sce)[,paste0("nmf",1:nFactors)])
seed1 = seed1>0
d1 = cbind.data.frame(cell.class=as.data.frame(colData(sce))[,"CellType.Final"],
                      seed1) %>% 
  group_by(cell.class) %>% add_tally(name="total") %>%
    group_by(cell.class, total) %>%
  summarise_at(paste0("nmf",1:nFactors), sum) %>%
  tidyr::pivot_longer(paste0("nmf",1:nFactors), values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)

seed2 = as.matrix(colData(sce)[,paste0("nmf",1:nFactors)])
seed2 = apply(seed2, 2, scale)
d2 = cbind.data.frame(cell.class=as.data.frame(colData(sce))[,"CellType.Final"],
                      seed2) %>% 
  group_by(cell.class) %>%
  summarise_at(paste0("nmf",1:nFactors), mean) %>% 
  tidyr::pivot_longer(paste0("nmf",1:nFactors), values_to="scaled.avg", names_to="nmf")

spf.ordered <- c("DRD1_MSN_C", "DRD1_MSN_A", "DRD2_MSN_A", "DRD2_MSN_B", "DRD1_MSN_B", "DRD1_MSN_D", 
  "Inh_A", "Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F",
  "Excitatory", "Astrocyte_A", "Astrocyte_B", "Oligo", "OPC", 
  "Microglia", "Endothelial", "Ependymal", "Neuron_Ambig")

dot.df = left_join(d1[,c("cell.class","nmf","prop")], 
                   d2[,c("cell.class","nmf","scaled.avg")]) %>%
  mutate(cell.class=factor(cell.class, 
                                     levels=spf.ordered))


nmf.ordered = c("nmf1", "nmf9", "nmf11", "nmf14", "nmf20", "nmf26", "nmf29", "nmf33", "nmf37", "nmf41", "nmf54", "nmf55", "nmf56", "nmf61", "nmf64", "nmf68", "nmf73", "nmf76", "nmf77", "nmf80", "nmf82", "nmf4", "nmf42", "nmf43", "nmf44", "nmf48", # General
                "nmf3", "nmf17", "nmf18", "nmf25", "nmf27", "nmf13", # General neuron
                "nmf19", "nmf21", "nmf23", # General MSN
                "nmf24", "nmf10", "nmf8", # DRD1 MSN C
                "nmf2", "nmf6", "nmf7", "nmf12", "nmf15", # DRD1/DRD2 MSN A
                "nmf5", # DRD2 MSN A
                "nmf39", # DRD2 MSN B
                "nmf34", "nmf36", "nmf31", # DRD1 MSN B
                "nmf46", # DRD1 MSN D
                "nmf32", # Inhib A
                "nmf66", # Inhib B
                "nmf50", # Inhib C
                "nmf72", # Inhib D
                "nmf63", # Inhib E
                "nmf74", #Inhib F
                "nmf65", "nmf67", # Excitatory
                "nmf47", "nmf70", "nmf53", "nmf22", "nmf35", "nmf62", # Astrocyte A/B
                "nmf28", "nmf38", "nmf40", "nmf45", "nmf49", "nmf51", "nmf52", "nmf58", "nmf59", "nmf60", #Oligo
                "nmf16", #OPC
                "nmf30", "nmf69", "nmf75","nmf79",  # Microglia
                "nmf78", "nmf81", # Endothelial
                "nmf57", # Ependymal
                "nmf71", "nmf83" # Neuron Ambig
                )
                
nmf.ordered.keep = setdiff(nmf.ordered, remove.nmf)
nmf.ordered.remove = intersect(nmf.ordered, remove.nmf)
sex_specific_factors <- rownames(corr_df)
nmf.ordered.keep <- nmf.ordered.keep[!nmf.ordered.keep %in% sex_specific_factors]
nmf.ordered.remove <- nmf.ordered.remove[!nmf.ordered.remove %in% sex_specific_factors]

dot.df$nmf_f = factor(dot.df$nmf, levels=c(rownames(corr_df), nmf.ordered.remove, nmf.ordered.keep))

pdf(file.path(plotDir, "NMF_by_cell_types.pdf"), width = 16, height = 6)
ggplot(dot.df, aes(x=nmf_f, y=cell.class, size=prop, color=scaled.avg))+
  geom_count()+theme_bw()+
  scale_size(range=c(0,3))+scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) + xlab("") + ylab("")
dev.off()

#spe dotplot
spe <- spe[ ,!is.na(spe[["domain"]])]
seed3 = as.matrix(colData(spe)[,nmf.ordered.keep])
seed3 = seed3>0
d3 = cbind.data.frame(domain=as.data.frame(colData(spe))[,"domain"],
                      seed3) %>% 
  group_by(domain) %>% add_tally(name="total") %>%
  group_by(domain, total) %>%
  summarise_at(nmf.ordered.keep, sum) %>%
  tidyr::pivot_longer(nmf.ordered.keep, values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)

seed4 = as.matrix(colData(spe)[,nmf.ordered.keep])
seed4 = apply(seed4, 2, scale)
d4 = cbind.data.frame(domain=as.data.frame(colData(spe))[,"domain"],
                      seed4) %>% 
  group_by(domain) %>%
  summarise_at(nmf.ordered.keep, mean) %>% 
  tidyr::pivot_longer(nmf.ordered.keep, values_to="scaled.avg", names_to="nmf")

dot.df2 = left_join(d3[,c("domain","nmf","prop")], 
                    d4[,c("domain","nmf","scaled.avg")])

domain.ordered <- c("MSN 1", "MSN 2", "MSN 3", "D1 islands", "Inhibitory", "Excitatory", 
  "WM", "Endothelial/Ependymal")

dot.df2$domain = factor(dot.df2$domain, levels = domain.ordered)
dot.df2$nmf_f = factor(dot.df2$nmf, levels=nmf.ordered.keep)

pdf(file.path(plotDir, "NMF_by_spatial_domains.pdf"), width = 12, height = 4)
ggplot(dot.df2, aes(x=nmf_f, y=domain, size=prop, color=scaled.avg))+
  geom_count()+theme_bw()+
  scale_size(range=c(0,3))+scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) + xlab("") + ylab("")
dev.off()