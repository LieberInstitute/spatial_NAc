# Load libraries
rm(list = ls())
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(scuttle)
library(edgeR)
library(dplyr)
library(here)

res_dir <- here("processed-data", "10_post_clustering_analysis", "02_spatial_registration_sn")
sn_registration <- readRDS(file = file.path(res_dir, "sn_cellType_registration.rds"))
enrichment_results <- sn_registration[["enrichment"]]
duplicated_genes <- enrichment_results[duplicated(enrichment_results$gene),]$gene

for (gene in duplicated_genes) {
    enrichment_results <- enrichment_results[-sample(which(enrichment_results$gene == gene), 1),]
}
aggregated <- data.frame(
    Astrocyte_A = enrichment_results$t_stat_Astrocyte_A,
    Astrocyte_B = enrichment_results$t_stat_Astrocyte_B,
    DRD1_MSN_A = enrichment_results$t_stat_DRD1_MSN_A,
    DRD1_MSN_B = enrichment_results$t_stat_DRD1_MSN_B,
    DRD1_MSN_C = enrichment_results$t_stat_DRD1_MSN_C,
    DRD1_MSN_D = enrichment_results$t_stat_DRD1_MSN_D,
    DRD2_MSN_A = enrichment_results$t_stat_DRD2_MSN_A, 
    DRD2_MSN_B   = enrichment_results$t_stat_DRD2_MSN_B, 
    Endothelial = enrichment_results$t_stat_Endothelial, 
    Ependymal = enrichment_results$t_stat_Ependymal, 
    Excitatory = enrichment_results$t_stat_Excitatory, 
    Inh_A = enrichment_results$t_stat_Inh_A, 
    Inh_B = enrichment_results$t_stat_Inh_B, 
    Inh_C = enrichment_results$t_stat_Inh_C, 
    Inh_D = enrichment_results$t_stat_Inh_D,
    Inh_E = enrichment_results$t_stat_Inh_E,
    Inh_F = enrichment_results$t_stat_Inh_F,
    Microglia = enrichment_results$t_stat_Microglia, 
    Oligo = enrichment_results$t_stat_Oligo, 
    OPC = enrichment_results$t_stat_OPC)

rownames(aggregated) <- enrichment_results$gene

write.table(aggregated, here::here("processed-data", "19_sLDSC_old","snRNA_seq","input_files", "snRNA_aggregated_de.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")
