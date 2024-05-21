#Completed in interactive job: srun --cpus-per-task=2 --mem=20G --pty --x11 bash
#cd cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc
#module load r_nac


library(SingleCellExperiment)
library(pheatmap)
library(reshape2)
library(here)

#Load NMF results and the sce object
nmf_res <- readRDS(here("processed-data","12_snRNA","NMF","NMF_Results.Rds"))
sce <-  readRDS(here("processed-data","12_snRNA","sce_CellType_noresiduals.Rds"))

sce
# class: SingleCellExperiment 
# dim: 36601 103785 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(103785): 1_AAACCCAAGACCAACG-1 1_AAACCCACAGTCAGCC-1 ...
# 20_TTTGTTGCAAGATGTA-1 20_TTTGTTGGTACGAAAT-1
# colData names(33): Sample Barcode ... sizeFactor CellType.Final
# reducedDimNames(5): GLMPCA_approx tSNE HARMONY tSNE_HARMONY
# umap_HARMONY
# mainExpName: NULL
# altExpNames(0):

###########Correlate Brain ID with NMF patterns#########
##onehot encoding of the brain variable##
#Make a dataframe of brain IDs 
data <- as.data.frame(sce$Brain_ID)
colnames(data) <- "Brain_ID"

head(data)
# Brain_ID
# 1   Br8325
# 2   Br8325
# 3   Br8325
# 4   Br8325
# 5   Br8325
# 6   Br8325

dim(data)
#[1] 103785      1

onehot_brain <-  dcast(data = data, rownames(data) ~ Brain_ID, length)

dim(onehot_brain)
#[1] 103785     11

head(onehot_brain)
# rownames(data) Br2720 Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325 Br8492
# 1              1      0      0      0      0      0      0      0      1      0
# 2             10      0      0      0      0      0      0      0      1      0
# 3            100      0      0      0      0      0      0      0      1      0
# 4           1000      0      0      0      0      0      0      0      1      0
# 5          10000      0      0      0      0      0      0      0      0      1
# 6         100000      0      0      0      0      0      0      0      0      0
# Br8667
# 1      0
# 2      0
# 3      0
# 4      0
# 5      0
# 6      1

#Make the rownames the first column
rownames(onehot_brain) <- onehot_brain[,1]

#Convert first column to numeric 
onehot_brain[,1]<-as.numeric(onehot_brain[,1])

#Reorder based on first column
onehot_brain <- onehot_brain[order(onehot_brain[,1],decreasing=FALSE),]

head(onehot_brain)
# rownames(data) Br2720 Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325 Br8492
# 1              1      0      0      0      0      0      0      0      1      0
# 2              2      0      0      0      0      0      0      0      1      0
# 3              3      0      0      0      0      0      0      0      1      0
# 4              4      0      0      0      0      0      0      0      1      0
# 5              5      0      0      0      0      0      0      0      1      0
# 6              6      0      0      0      0      0      0      0      1      0
# Br8667
# 1      0
# 2      0
# 3      0
# 4      0
# 5      0
# 6      0

#Remove the first column
onehot_brain[,1] <- NULL

#Correlate brain ID with NMF patterns. 
pdf(here("plots","12_snRNA","NMF", "NMF_BrainID_correlation_heatmap.pdf"))
pheatmap(cor(t(nmf_res@h),onehot_brain), fontsize_row = 5)
dev.off()

###########Correlate QC measures with NMF patterns#########
##onehot encoding of the brain variable##
#Make a dataframe of brain IDs 
QC_vars <- colData(sce)[,c("detected","sum","subsets_Mito_percent")]

head(QC_vars)
# DataFrame with 6 rows and 3 columns
# detected       sum subsets_Mito_percent
# <integer> <numeric>            <numeric>
# 1_AAACCCAAGACCAACG-1      1744      3673            0.0544514
# 1_AAACCCACAGTCAGCC-1      2440      5977            0.0000000
# 1_AAACCCATCACCCTCA-1      6355     21932            0.0775123
# 1_AAACCCATCTACGGGC-1      1845      3174            0.1260239
# 1_AAACGAAAGGTGAGCT-1      1918      3560            0.0280899
1_AAACGAACATAGACTC-1      6795     27502            0.0145444

#Make the detected numeric
QC_vars$detected <- as.numeric(QC_vars$detected)

#Convert to matrix
QC_vars <- as.matrix(QC_vars)

pdf(here("plots","12_snRNA","NMF","NMF_QC_correlation_heatmap.pdf"))
pheatmap(cor(t(nmf_res@h),QC_vars), fontsize_row = 5)
dev.off()


###########Correlate CellType with NMF patterns#########
#Create datafarme of celltype 
ct_data <- as.data.frame(sce$CellType.Final)
colnames(ct_data) <- "CellType.Final" 

#One hot encode the cell type
onehot_CellType <-  dcast(data = ct_data, rownames(ct_data) ~ CellType.Final, length)

head(onehot_CellType)
# rownames(ct_data) Astrocyte_A Astrocyte_B DRD1_MSN_A DRD1_MSN_B DRD1_MSN_C
# 1                 1           0           0          0          0          0
# 2                10           0           0          0          0          1
# 3               100           0           0          1          0          0
# 4              1000           0           0          0          0          0
# 5             10000           0           0          0          0          0
# 6            100000           0           0          0          0          0
# DRD1_MSN_D DRD2_MSN_A DRD2_MSN_B Ependymal Excitatory Inh_A Inh_B Inh_C Inh_D
# 1          0          0          0         0          0     0     0     0     0
# 2          0          0          0         0          0     0     0     0     0
# 3          0          0          0         0          0     0     0     0     0
# 4          0          0          1         0          0     0     0     0     0
# 5          0          0          1         0          0     0     0     0     0
# 6          0          0          1         0          0     0     0     0     0
# Inh_E Inh_F Inh_G Inh_H Inh_I Macrophage Microglia
# 1     0     0     0     0     0          0         0
# 2     0     0     0     0     0          0         0
# 3     0     0     0     0     0          0         0
# 4     0     0     0     0     0          0         0
# 5     0     0     0     0     0          0         0
# 6     0     0     0     0     0          0         0
# Mural_Endothelial_Fibroblast Oligo OPC T-Cell
# 1                            0     1   0      0
# 2                            0     0   0      0
# 3                            0     0   0      0
# 4                            0     0   0      0
# 5                            0     0   0      0
# 6                            0     0   0      0





