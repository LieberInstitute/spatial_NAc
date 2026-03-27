library(SpatialExperiment)
library(escheR)
library(Seurat)
library(here)

spe <- readRDS(here("processed-data","xenium","SPEs","spe_NormCounts.Rds"))

dim(spe)

spe

#load in the snRNA-seq predictions
predictions <- readRDS(here("processed-data","12_snRNA","xenium_snRNA_predictions.Rds"))

dim(predictions)

#Sanity checks before adding 
stopifnot(identical(nrow(predictions),ncol(spe)))


stopifnot(identical(rownames(predictions),rownames(colData(spe))))


stopifnot(identical(rownames(predictions),colnames(spe)))


#Prediction dataframes have the predicted id and the weighted probability score for each designation. 
#Add that information to the colData of the spe object. prediction.score.max
colData(spe)$snRNA_predicted_CellType <- predictions[,"predicted.id"]
colData(spe)$predicted_max_score <- predictions[,"prediction.score.max"]

table(colData(spe)$snRNA_predicted_CellType)


#load cluster colors from snRNA-seq paper. 
load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda",
     verbose = TRUE)

#Remove neuron ambig
cluster_cols <- cluster_cols[-14]
cluster_cols


spe$snRNA_predicted_CellType <- factor(x = spe$snRNA_predicted_CellType,
                                       levels = c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                                                  "DRD2_MSN_A","DRD2_MSN_B",
                                                  "Inh_A","Inh_B","Inh_C","Inh_D","Inh_E","Inh_F",
                                                  "Excitatory",
                                                  "Oligo","OPC",
                                                  "Astrocyte_A","Astrocyte_B","Ependymal",
                                                  "Microglia",
                                                  "Endothelial"))


library(ggplot2)

max_score_box <- ggplot(colData(spe),aes(x = snRNA_predicted_CellType, 
                                         y = predicted_max_score,
                                         fill = snRNA_predicted_CellType)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.50) +
  theme_bw() +
  scale_fill_manual(values = cluster_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "Predicted Cell Type",
       y = "Max Score")
ggsave(max_score_box,filename = here("plots","xenium","Label_Transfer",
                                     "Max_score_boxplot.pdf"))

#Plot the predicted IDs on the tissue sections .
for(i in unique(spe$Sample)){
  print(i)
  spe_sub <- spe[,spe$Sample == i]
  p <- make_escheR(spe_sub) |>
    add_fill("snRNA_predicted_CellType") +
    scale_fill_manual(values = cluster_cols)
  ggsave(filename = here("plots","xenium","Label_Transfer",paste0(i,".png")),
         height = 25, width = 25)
}

#Plot the predicted IDs on the tissue sections .
for(sample in unique(spe$Sample)){
  print(sample)
  spe_sub <- spe[,spe$Sample == sample]
  dir.create(path = here("plots","xenium","Label_Transfer","Cluster_Specific",sample))
  for(cluster in unique(spe$snRNA_predicted_CellType)){
    print(cluster)
    spe_sub$Cluster <- ifelse(spe_sub$snRNA_predicted_CellType == cluster,
                          cluster,
                          "Other")
    p <- make_escheR(spe_sub) |>
      add_fill("Cluster") +
      scale_fill_manual(values = c(cluster_cols[cluster], "Other" = "grey80"))
    ggsave(filename = here("plots","xenium","Label_Transfer","Cluster_Specific",sample,paste0(cluster,"_",sample,".png")),
           height = 25, width = 25)
  }
}


#Plot some important genes. 
genes <- c("ADORA2A","CALB1","PPP1R1B","CARTPT","ARHGAP36","PDYN","PENK","OPRK1","DRD1","DRD2","FOXJ1","GABRQ","RXFP1")

#Plot genes to help identify WM tracts and NAc proper
spe$ADORA2A <- assay(spe,"nucleus_normcounts")["ADORA2A",]
spe$CALB1 <- assay(spe,"nucleus_normcounts")["CALB1",]
spe$PPP1R1B <- assay(spe,"nucleus_normcounts")["PPP1R1B",]
spe$CARTPT <- assay(spe,"nucleus_normcounts")["CARTPT",]
spe$ARHGAP36 <- assay(spe,"nucleus_normcounts")["ARHGAP36",]
spe$PDYN <- assay(spe,"nucleus_normcounts")["PDYN",]
spe$PENK <- assay(spe,"nucleus_normcounts")["PENK",]
spe$OPRK1 <- assay(spe,"nucleus_normcounts")["OPRK1",]
spe$DRD1 <- assay(spe,"nucleus_normcounts")["DRD1",]
spe$DRD2 <- assay(spe,"nucleus_normcounts")["DRD2",]
spe$FOXJ1 <- assay(spe,"nucleus_normcounts")["FOXJ1",]
spe$GABRQ <- assay(spe,"nucleus_normcounts")["GABRQ",]
spe$RXFP1 <- assay(spe,"nucleus_normcounts")["RXFP1",]



for(sample in unique(spe$Sample)){
    print(sample)
    sub_spe <- spe[,spe$Sample == sample]
    for(gene in genes){
        print(gene)
        p1 <- make_escheR(sub_spe) |>
                add_fill(gene) +
                scale_fill_gradientn(colors = c("grey80","orange","red"))
        ggsave(plot = p1,
               filename = here("plots","xenium","Expression",
                              paste0(sample,"_",gene,".png")),
               height = 18,width = 18)
  }
}



###Reproduciblity
sessionInfo()
