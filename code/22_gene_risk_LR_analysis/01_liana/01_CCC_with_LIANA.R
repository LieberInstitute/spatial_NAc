library(tidyverse)
library(magrittr)
library(liana)
library(here)
library(SingleCellExperiment)

# Read in the single cell data
dat_dir <- here("processed-data", "12_snRNA")
sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))

# Remove neuron ambig, and add the cell type labels as sce idents
sce <- sce[ ,!sce$CellType.Final == "Neuron_Ambig"]
colData(sce)$CellType.Final <- factor(colData(sce)$CellType.Final)
colLabels(sce) <- colData(sce)$CellType.Final

# Change the rownames to gene names because it's consistent with Omnipath
rownames(sce) <- rowData(sce)$gene_name

print("The dimensions of the SCE object are:")
print(dim(sce))

# Run liana
liana_results <- liana_wrap(sce)
res_dir <- here("processed-data", "22_gene_risk_LR_analysis", "01_liana")

saveRDS(liana_results, file = file.path(res_dir, "liana_results.rds"))

# Consensus ranks
liana_consensus <- liana_results %>% liana_aggregate()
saveRDS(liana_consensus,
 file = file.path(res_dir, "liana_consensus.rds"))

sessionInfo()