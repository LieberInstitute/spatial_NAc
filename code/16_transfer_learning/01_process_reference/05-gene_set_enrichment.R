####code for obtaining gene set enrichment for select NMF patterns
library(dplyr)
library(ggplot2)
library(fgsea)
library(reactome.db)
library(org.Hs.eg.db)

set.seed(123)

# Load NMFs
#spec <- matrix(
#    c(  "data", "d", 1, "character", "Specify the dataset to be used?"
#    ),
#    byrow = TRUE, ncol = 5
#)
#opt <- getopt(spec)

opt <- list()
opt$data <- "human_NAc"
print(opt$data)

# Read data and create Seurat object
dat_dir <- here::here("processed-data", "12_snRNA")
res_dir <- here::here("processed-data", "16_transfer_learning","01_process_reference", "RCppML", opt$data)
plot_dir <- here::here("plots", "16_transfer_learning","01_process_reference", "RCppML", opt$data)

#if(opt$data == "human_NAc"){
#  sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
#}else{
#  if(opt$data == "rat_case_control"){
#    sce <- readRDS(file = file.path(dat_dir, "NAc_Combo_Integrated.RDS"))
#  }else{
#        stop("Invalid input data set")
#  }
#}

x <- readRDS(file = file.path(res_dir,paste0("nmf_results.rds")))
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

# Prep reactome
xx <- as.list(reactomePATHID2NAME)
reactome.h <- xx[grep("^Homo",xx)]
x <- as.list(reactomePATHID2EXTID)
reactome.h = reactome.h[intersect(names(reactome.h), names(x))]
x.h <- x[names(reactome.h)]
identical(names(x.h), names(reactome.h))
reactome.h = gsub("Homo sapiens: ","",reactome.h)
names(x.h) = reactome.h

################# GSEA
non0.nmf3 = rownames(loads)[loads[,"nmf3"]>0]
non0.3.id = mapIds(org.Hs.eg.db, keys=non0.nmf3, keytype="ENSEMBL", column="ENTREZID", multiVals = "first")
names(non0.3.id) = non0.nmf3
non0.3.id = non0.3.id[!is.na(non0.3.id)]
pathways.3 <- reactomePathways(non0.3.id)
pathways.3 <- x.h[names(pathways.3)]

nmf3.stats = loads[names(non0.3.id),"nmf3"]
names(nmf3.stats) = non0.3.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
nmf.results.3 = fgseaMultilevel(pathways.3, stats=nmf3.stats, scoreType="pos", minSize=15, maxSize=500)
nmf.results.3$leadingEdge2 = sapply(nmf.results.3$leadingEdge, paste, collapse="/")
write.csv(nmf.results.3[,c(1:7,9)], file.path(res_dir ,"nmf3_reactome_results.csv"))

# Repeat for NMF 4
non0.nmf4 = rownames(loads)[loads[,"nmf4"]>0]
non0.4.id = mapIds(org.Hs.eg.db, keys=non0.nmf4, keytype="ENSEMBL", column="ENTREZID", multiVals = "first")
names(non0.4.id) = non0.nmf4
non0.4.id = non0.4.id[!is.na(non0.4.id)]
pathways.4 <- reactomePathways(non0.4.id)
pathways.4 <- x.h[names(pathways.4)]

nmf4.stats = loads[names(non0.4.id),"nmf4"]
names(nmf4.stats) = non0.4.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
nmf.results.4 = fgseaMultilevel(pathways.4, stats=nmf4.stats, scoreType="pos", minSize=15, maxSize=500)
nmf.results.4$leadingEdge2 = sapply(nmf.results.4$leadingEdge, paste, collapse="/")
write.csv(nmf.results.4[,c(1:7,9)], file.path(res_dir ,"nmf4_reactome_results.csv"))

# Repeat for NMF 7
non0.nmf7 = rownames(loads)[loads[,"nmf7"]>0]
non0.7.id = mapIds(org.Hs.eg.db, keys=non0.nmf7, keytype="ENSEMBL", column="ENTREZID", multiVals = "first")
names(non0.7.id) = non0.nmf7
non0.7.id = non0.7.id[!is.na(non0.7.id)]
pathways.7 <- reactomePathways(non0.7.id)
pathways.7 <- x.h[names(pathways.7)]

nmf7.stats = loads[names(non0.7.id),"nmf7"]
names(nmf7.stats) = non0.7.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
nmf.results.7 = fgseaMultilevel(pathways.7, stats=nmf7.stats, scoreType="pos", minSize=15, maxSize=500)
nmf.results.7$leadingEdge2 = sapply(nmf.results.7$leadingEdge, paste, collapse="/")
write.csv(nmf.results.7[,c(1:7,9)], file.path(res_dir ,"nmf7_reactome_results.csv"))

# Repeat for NMF 10
non0.nmf10 = rownames(loads)[loads[,"nmf10"]>0]
non0.10.id = mapIds(org.Hs.eg.db, keys=non0.nmf10, keytype="ENSEMBL", column="ENTREZID", multiVals = "first")
names(non0.10.id) = non0.nmf10
non0.10.id = non0.10.id[!is.na(non0.10.id)]
pathways.10 <- reactomePathways(non0.10.id)
pathways.10 <- x.h[names(pathways.10)]

nmf10.stats = loads[names(non0.10.id),"nmf10"]
names(nmf10.stats) = non0.10.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
nmf.results.10 = fgseaMultilevel(pathways.10, stats=nmf10.stats, scoreType="pos", minSize=15, maxSize=500)
nmf.results.10$leadingEdge2 = sapply(nmf.results.10$leadingEdge, paste, collapse="/")
write.csv(nmf.results.10[,c(1:7,9)], file.path(res_dir ,"nmf10_reactome_results.csv"))

# Repeat for NMF 34
non0.nmf34 = rownames(loads)[loads[,"nmf34"]>0]
non0.34.id = mapIds(org.Hs.eg.db, keys=non0.nmf34, keytype="ENSEMBL", column="ENTREZID", multiVals = "first")
names(non0.34.id) = non0.nmf34
non0.34.id = non0.34.id[!is.na(non0.34.id)]
pathways.34 <- reactomePathways(non0.34.id)
pathways.34 <- x.h[names(pathways.34)]

nmf34.stats = loads[names(non0.34.id),"nmf34"]
names(nmf34.stats) = non0.34.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
nmf.results.34 = fgseaMultilevel(pathways.34, stats=nmf34.stats, scoreType="pos", minSize=15, maxSize=500)
nmf.results.34$leadingEdge2 = sapply(nmf.results.34$leadingEdge, paste, collapse="/")
write.csv(nmf.results.34[,c(1:7,9)], file.path(res_dir ,"nmf34_reactome_results.csv"))

# Repeat for NMF 35
non0.nmf35 = rownames(loads)[loads[,"nmf35"]>0]
non0.35.id = mapIds(org.Hs.eg.db, keys=non0.nmf35, keytype="ENSEMBL", column="ENTREZID", multiVals = "first")
names(non0.35.id) = non0.nmf35
non0.35.id = non0.35.id[!is.na(non0.35.id)]
pathways.35 <- reactomePathways(non0.35.id)
pathways.35 <- x.h[names(pathways.35)]

nmf35.stats = loads[names(non0.35.id),"nmf35"]
names(nmf35.stats) = non0.35.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
nmf.results.35 = fgseaMultilevel(pathways.35, stats=nmf35.stats, scoreType="pos", minSize=15, maxSize=500)
nmf.results.35$leadingEdge2 = sapply(nmf.results.35$leadingEdge, paste, collapse="/")
write.csv(nmf.results.35[,c(1:7,9)], file.path(res_dir ,"nmf35_reactome_results.csv"))

# Repeat for NMF 44
non0.nmf44 = rownames(loads)[loads[,"nmf44"]>0]
non0.44.id = mapIds(org.Hs.eg.db, keys=non0.nmf44, keytype="ENSEMBL", column="ENTREZID", multiVals = "first")
names(non0.44.id) = non0.nmf44
non0.44.id = non0.44.id[!is.na(non0.44.id)]
pathways.44 <- reactomePathways(non0.44.id)
pathways.44 <- x.h[names(pathways.44)]

nmf44.stats = loads[names(non0.44.id),"nmf44"]
names(nmf44.stats) = non0.44.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
nmf.results.44 = fgseaMultilevel(pathways.44, stats=nmf44.stats, scoreType="pos", minSize=15, maxSize=500)
nmf.results.44$leadingEdge2 = sapply(nmf.results.44$leadingEdge, paste, collapse="/")
write.csv(nmf.results.44[,c(1:7,9)], file.path(res_dir ,"nmf44_reactome_results.csv"))


# GSEA tables

nmf10.terms = c("Signaling by GPCR","G alpha (i) signalling events","G alpha (q) signalling events", "Opioid Signalling")
pdf(file.path(plot_dir, "GSEA_table_nmf10.pdf"), width = 7, height = 3)
plotGseaTable(x.h[nmf10.terms],
              nmf10.stats, nmf.results.10, 
              gseaParam = 0.5)
dev.off()

nmf7.terms = c("Signaling by GPCR", "Opioid Signalling", "Activation of NMDA receptors and postsynaptic events")
pdf(file.path(plot_dir, "GSEA_table_nmf7.pdf"), width = 7, height = 3)
plotGseaTable(x.h[nmf7.terms],
              nmf7.stats, nmf.results.7, 
              gseaParam = 0.5)
dev.off()

nmf4.terms = c("Signaling by GPCR", "G alpha (i) signalling events", "Hedgehog ligand biogenesis")
pdf(file.path(plot_dir, "GSEA_table_nmf4.pdf"), width = 7, height = 3)
plotGseaTable(x.h[nmf4.terms],
              nmf4.stats, nmf.results.4, 
              gseaParam = 0.5)
dev.off()

nmf3.terms = c("Signaling by GPCR", "Opioid Signalling", "Sensory Perception", "Neurexins and neuroligins")
pdf(file.path(plot_dir, "GSEA_table_nmf3.pdf"), width = 7, height = 3)
plotGseaTable(x.h[nmf3.terms],
              nmf3.stats, nmf.results.3, 
              gseaParam = 0.5)
dev.off()

nmf34.terms = c("Activation of NMDA receptors and postsynaptic events", 
"Activation of kainate receptors upon glutamate binding", "DAG and IP3 signaling")
pdf(file.path(plot_dir, "GSEA_table_nmf34.pdf"), width = 7, height = 3)
plotGseaTable(x.h[nmf34.terms],
              nmf34.stats, nmf.results.34, 
              gseaParam = 0.5)
dev.off()

nmf35.terms = c("Activation of NMDA receptors and postsynaptic events", 
"Axon guidance", "Glutamate binding, activation of AMPA receptors and synaptic plasticity", 
"Regulation of insulin secretion")
pdf(file.path(plot_dir, "GSEA_table_nmf35.pdf"), width = 7, height = 3)
plotGseaTable(x.h[nmf35.terms],
              nmf35.stats, nmf.results.35, 
              gseaParam = 0.5)
dev.off()

nmf44.terms = c("Activation of NMDA receptors and postsynaptic events", 
"Cell junction organization", "Dopamine Neurotransmitter Release Cycle", "Interaction between L1 and Ankyrins")
pdf(file.path(plot_dir, "GSEA_table_nmf44.pdf"), width = 7, height = 3)
plotGseaTable(x.h[nmf44.terms],
              nmf44.stats, nmf.results.44, 
              gseaParam = 0.5)
dev.off()

