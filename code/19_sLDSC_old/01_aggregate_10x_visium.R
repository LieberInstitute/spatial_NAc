# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(scuttle)
library(edgeR)
library(dplyr)
library(here)

# Read in the modeling results
modeling_rdata_path = here(
    'processed-data', '10_post_clustering_analysis','01_pseudobulk_markers', '01_precast', 'pseudobulk_capture_area', 'final_clusters', 'model_results_precast_clusters.Rdata')
load(modeling_rdata_path)

enrichment_results <- modeling_results[["enrichment"]]
duplicated_genes <- enrichment_results[duplicated(enrichment_results$gene),]$gene

for (gene in duplicated_genes) {
    enrichment_results <- enrichment_results[-sample(which(enrichment_results$gene == gene), 1),]
}

aggregated <- data.frame(
    D1_islands = enrichment_results$t_stat_D1.islands,
    Endothelial_Ependymal = enrichment_results$t_stat_Endothelial.Ependymal,
    Excitatory = enrichment_results$t_stat_Excitatory,
    Inhibitory = enrichment_results$t_stat_Inhibitory,
    MSN_1 = enrichment_results$t_stat_MSN.1,
    MSN_2 = enrichment_results$t_stat_MSN.2,
    MSN_3 = enrichment_results$t_stat_MSN.3, 
    WM   = enrichment_results$t_stat_WM)

rownames(aggregated) <- enrichment_results$gene

write.table(aggregated, here::here("processed-data", "19_sLDSC_old","10x_visium","input_files", "visium_aggregated_de.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")
