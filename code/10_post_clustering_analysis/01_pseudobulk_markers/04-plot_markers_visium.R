library(here)
library(PRECAST)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(spatialNAcUtils)
library(purrr)
library(ggpubr)
library(ggsci)
library(dittoSeq)
library(getopt)
library(pheatmap)
library(cowplot)

plotDir <- here("plots", "10_post_clustering_analysis", "01_pseudobulk_markers", 
                "01_precast", "plot_select_markers")
subregion_genes <- list(
    "shell" = c(
        "CALB2", "CPNE4", "ARHGAP36", "OPRK1", "HTR4", "CARTPT", "PYDC1", "PNMA5", "CALCR", "CBLN4"
    ),
    "core" = toupper(
        c(
            "CALB1", "LAMP5", "PEG10", "PDE1A", "ADORA2A", "PTPN7", "PDE10A"
        )
    ),
    "D1_islands" = c("FOXP2", "OPRM1", "RXFP1", "CPNE4", "DRD1", "PDYN", "PENK"), 
    "white_matter" = c("MBP", "MOBP", "PPP1R1B", "OLIG1"), 
    "inhibitory" = c("KIT", "SST", "NPY", "CHODL")
)

select_samples <- c("Br6432", "Br6522", "Br3942")

spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5"
)
spe <- loadHDF5SummarizedExperiment(spe_dir)
stopifnot(all(unlist(subregion_genes) %in% rowData(spe)$gene_name))

# Subset to select samples
spe <- spe[ ,spe$donor %in% select_samples]
spe$exclude_overlapping[spe$sample_id_original == "V11D01-384_A1" & spe$overlap_slide == "V11D01-384_D1"] <- FALSE
spe$exclude_overlapping[spe$sample_id_original == "V11D01-384_D1" & spe$overlap_slide == "V11D01-384_A1"] <- TRUE
# Read in clustering results
precast_path = here(
    'processed-data', '07_spatial_domains', '01_precast', 'nnSVG_precast', 'final_clusters', 'precast_clusters.csv')

cluster_col <-"spatial_domains"
spe[[cluster_col]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(precast_path), by = 'key') |>
    pull(cluster) |>
    as.factor()

spe <- spe[ ,!is.na(spe[[cluster_col]])]
rownames(spe) <- rowData(spe)$gene_name


for(i in c(1:length(subregion_genes$D1_islands))){
    cat(subregion_genes$D1_islands[i], "\n")
    p_D1islands <- list()
    for(donor in select_samples){
        p_D1islands[[donor]] <- vis_gene(spe, sampleid = donor, geneid = subregion_genes$D1_islands[i], is_stitched = TRUE)
    }
    dir.create(paste0(plotDir, "/D1_islands"))
    pdf(file.path(plotDir, "D1_islands", paste0(subregion_genes$D1_islands[i], ".pdf")), width = 6, height = 6)
    p_D1islands[["Br6432"]]
    p_D1islands[["Br6522"]]
    p_D1islands[["Br3942"]]
    dev.off()
}

dir.create(paste0(plotDir, "/core"))
for(i in c(1:length(subregion_genes$core))){
    cat(subregion_genes$core[i], "\n")
    p_core <- list()
    for(donor in select_samples){
        p_core[[donor]] <- vis_gene(spe, sampleid = donor, geneid = subregion_genes$core[i], is_stitched = TRUE)
    }
    pdf(file.path(plotDir, "core", paste0(subregion_genes$core[i], ".pdf")), width = 6, height = 6)
    p_core[["Br6432"]]
    p_core[["Br6522"]]
    p_core[["Br3942"]]
    dev.off()
}

dir.create(paste0(plotDir, "/shell"))
for(i in c(1:length(subregion_genes$shell))){
    cat(subregion_genes$shell[i], "\n")
    p_shell <- list()
    for(donor in select_samples){
        p_shell[[donor]] <- vis_gene(spe, sampleid = donor, geneid = subregion_genes$shell[i], is_stitched = TRUE)
    }
    pdf(file.path(plotDir, "shell", paste0(subregion_genes$shell[i], ".pdf")), width = 6, height = 6)
    p_shell[["Br6432"]]
    p_shell[["Br6522"]]
    p_shell[["Br3942"]]
    dev.off()
}
dir.create(paste0(plotDir, "/white_matter"))
for(i in c(1:length(subregion_genes$white_matter))){
    p_white_matter <- list()
    for(donor in select_samples){
        p_white_matter[[donor]] <- vis_gene(spe, sampleid = donor, geneid = subregion_genes$white_matter[i], is_stitched = TRUE)
    }
    pdf(file.path(plotDir, "white_matter", paste0(subregion_genes$white_matter[i], ".pdf")), width = 6, height = 6)
    p_white_matter[["Br6432"]]
    p_white_matter[["Br6522"]]
    p_white_matter[["Br3942"]]
    dev.off()
}

dir.create(paste0(plotDir, "/inhibitory"))
for(i in c(1:length(subregion_genes$inhibitory))){
    p_inhibitory <- list()
    for(donor in select_samples){
        p_inhibitory[[donor]] <- vis_gene(spe, sampleid = donor, geneid = subregion_genes$inhibitory[i], is_stitched = TRUE)
    }
    pdf(file.path(plotDir, "inhibitory", paste0(subregion_genes$inhibitory[i], ".pdf")), width = 6, height = 6)
    p_inhibitory[["Br6432"]]
    p_inhibitory[["Br6522"]]
    p_inhibitory[["Br3942"]]
    dev.off()
}

pdf(file.path(plotDir, "ADORA2A.pdf"), width = 6, height = 6)
vis_gene(spe, sampleid = "Br6522", geneid = "ADORA2A", is_stitched = TRUE) + ggtitle("ADORA2A") + theme(plot.title = element_text(face = "bold.italic"), legend.title = element_text(size = 14))
dev.off()

pdf(file.path(plotDir, "ARHGAP36.pdf"), width = 6, height = 6)
vis_gene(spe, sampleid = "Br6522", geneid = "ARHGAP36", is_stitched = TRUE) + ggtitle("ARHGAP36") + theme(plot.title = element_text(face = "bold.italic"), legend.title = element_text(size = 14))
dev.off()

pdf(file.path(plotDir, "PDE10A.pdf"), width = 6, height = 6)
vis_gene(spe, sampleid = "Br6522", geneid = "PDE10A", is_stitched = TRUE) + ggtitle("PDE10A") + theme(plot.title = element_text(face = "bold.italic"), legend.title = element_text(size = 14))
dev.off()

pdf(file.path(plotDir, "CNR1.pdf"), width = 6, height = 6)
vis_gene(spe, sampleid = "Br6522", geneid = "CNR1", is_stitched = TRUE) + ggtitle("CNR1") + theme(plot.title = element_text(face = "bold.italic"), legend.title = element_text(size = 14))
dev.off()

pdf(file.path(plotDir, "NPY.pdf"), width = 6, height = 6)
vis_gene(spe, sampleid = "Br6522", geneid = "NPY", is_stitched = TRUE) + ggtitle("NPY") + theme(plot.title = element_text(face = "bold.italic"), legend.title = element_text(size = 14))
dev.off()

pdf(file.path(plotDir, "CALB1.pdf"), width = 6, height = 6)
vis_gene(spe, sampleid = "Br6522", geneid = "CALB1", is_stitched = TRUE) + ggtitle("CALB1") + theme(plot.title = element_text(face = "bold.italic"), legend.title = element_text(size = 14))
dev.off()

pdf(file.path(plotDir, "CARTPT.pdf"), width = 6, height = 6)
vis_gene(spe, sampleid = "Br6522", geneid = "CARTPT", is_stitched = TRUE) + ggtitle("CARTPT") + theme(plot.title = element_text(face = "bold.italic"), legend.title = element_text(size = 14))
dev.off()

pdf(file.path(plotDir, "MBP.pdf"), width = 6, height = 6)
vis_gene(spe, sampleid = "Br6522", geneid = "MBP", is_stitched = TRUE) + ggtitle("MBP") + theme(plot.title = element_text(face = "bold.italic"), legend.title = element_text(size = 14))
dev.off()

pdf(file.path(plotDir, "PPP1R1B.pdf"), width = 6, height = 6)
vis_gene(spe, sampleid = "Br6522", geneid = "PPP1R1B", is_stitched = TRUE) + ggtitle("PPP1R1B") + theme(plot.title = element_text(face = "bold.italic"), legend.title = element_text(size = 14))
dev.off()

pdf(file.path(plotDir, "DRD1.pdf"), width = 6, height = 6)
vis_gene(spe, sampleid = "Br6522", geneid = "DRD1", is_stitched = TRUE) + ggtitle("DRD1") + theme(plot.title = element_text(face = "bold.italic"), legend.title = element_text(size = 14))
dev.off()

pdf(file.path(plotDir, "DRD2.pdf"), width = 6, height = 6)
vis_gene(spe, sampleid = "Br6522", geneid = "DRD2", is_stitched = TRUE) + ggtitle("DRD2") + theme(plot.title = element_text(face = "bold.italic"), legend.title = element_text(size = 14))
dev.off()

pdf(file.path(plotDir, "GABRQ.pdf"), width = 6, height = 6)
vis_gene(spe, sampleid = "Br6522", geneid = "GABRQ", is_stitched = TRUE) + ggtitle("GABRQ") + theme(plot.title = element_text(face = "bold.italic"), legend.title = element_text(size = 14))
dev.off()

pdf(file.path(plotDir, "RXFP1.pdf"), width = 6, height = 6)
vis_gene(spe, sampleid = "Br6522", geneid = "RXFP1", is_stitched = TRUE) + ggtitle("RXFP1") + theme(plot.title = element_text(face = "bold.italic"), legend.title = element_text(size = 14))
dev.off()

pdf(file.path(plotDir, "OPRM1.pdf"), width = 6, height = 6)
vis_gene(spe, sampleid = "Br6522", geneid = "OPRM1", is_stitched = TRUE) + ggtitle("OPRM1") + theme(plot.title = element_text(face = "bold.italic"), legend.title = element_text(size = 14))
dev.off()

pdf(file.path(plotDir, "FOXP2.pdf"), width = 6, height = 6)
vis_gene(spe, sampleid = "Br6522", geneid = "FOXP2", is_stitched = TRUE) + ggtitle("FOXP2") + theme(plot.title = element_text(face = "bold.italic"), legend.title = element_text(size = 14))
dev.off()