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

spot_plot2 <- function(spe, sample_id, image_id = "lowres",
    title = sprintf("%s_%s", sample_id, var_name), var_name,
    multi_gene_method = c("z_score", "pca", "sparsity"),
    include_legend = TRUE, is_discrete, colors = NULL,
    assayname = "logcounts", minCount = 0.5, spatial = FALSE) {
    #   This value was determined empirically, and results in good spot sizes.
    #   Note that it's sample-independent, and the final spot size to pass to
    #   'vis_gene' or 'vis_clus' uses this value along with the image
    #   dimensions, scale factors, and spatial coordinates for this particular
    #   sample
    IDEAL_POINT_SIZE <- 200

    ############################################################################
    #   Check validity of arguments
    ############################################################################

    # (Note that 'sample_id', 'var_name', 'assayname', 'minCount', and 'colors'
    # are not checked for validity here, since spatialLIBD functions handle
    # their validity)

    multi_gene_method <- rlang::arg_match(multi_gene_method)

    #   Check validity of spatial coordinates
    if (!all(c("pxl_col_in_fullres", "pxl_row_in_fullres") == sort(colnames(spatialCoords(spe))))) {
        stop("Abnormal spatial coordinates: should have 'pxl_row_in_fullres' and 'pxl_col_in_fullres' columns.")
    }

    #   State assumptions about columns expected to be in the colData
    expected_cols <- c("array_row", "array_col", "sample_id", "exclude_overlapping")
    if (!all(expected_cols %in% colnames(colData(spe)))) {
        stop(
            sprintf(
                'Missing at least one of the following colData columns: "%s"',
                paste(expected_cols, collapse = '", "')
            )
        )
    }

    #   Subset to specific sample ID, and ensure overlapping spots are dropped
    subset_cols <- (spe$sample_id == sample_id) &
        (is.na(spe$exclude_overlapping) | !spe$exclude_overlapping)
    if (length(which(subset_cols)) == 0) {
        stop("No non-excluded spots belong to this sample. Perhaps check spe$exclude_overlapping for issues.")
    }
    spe_small <- spe[, subset_cols]

    ############################################################################
    #   Compute an appropriate spot size for this sample
    ############################################################################

    #   Determine some pixel values for the horizontal bounds of the spots
    MIN_COL <- min(spatialCoords(spe_small)[, "pxl_row_in_fullres"])
    MAX_COL <- max(spatialCoords(spe_small)[, "pxl_row_in_fullres"])

    #   The distance between spots (in pixels) is double the average distance
    #   between array columns
    INTER_SPOT_DIST_PX <- 2 * (MAX_COL - MIN_COL) /
        (max(spe_small$array_col) - min(spe_small$array_col))

    #   Find the appropriate spot size for this donor. This can vary because
    #   ggplot downscales a plot the fit desired output dimensions (in this
    #   case presumably a square region on a PDF), and stitched images can vary
    #   in aspect ratio. Also, lowres images always have a larger image
    #   dimension of 1200, no matter how many spots fit in either dimension.
    small_image_data <- imgData(spe_small)[
        imgData(spe_small)$image_id == image_id,
    ]
    spot_size <- IDEAL_POINT_SIZE * INTER_SPOT_DIST_PX *
        small_image_data$scaleFactor / max(dim(small_image_data$data[[1]]))


    ############################################################################
    #   Produce the plot
    ############################################################################

    #   If the quantity to plot is discrete, use 'vis_clus'. Otherwise use
    #   'vis_gene'.
    if (is_discrete) {
        #   Supply a color scale if 'color' is not NULL. Otherwise, fall back
        #   upon 'vis_clus' defaults
        if (is.null(colors)) {
            p <- vis_clus(
                spe_small,
                sampleid = sample_id, image_id = image_id, clustervar = var_name, auto_crop = FALSE,
                return_plots = TRUE, spatial = spatial, point_size = spot_size
            )
        } else {
            p <- vis_clus(
                spe_small,
                sampleid = sample_id, image_id = image_id, clustervar = var_name, auto_crop = FALSE,
                return_plots = TRUE, spatial = spatial, colors = colors,
                point_size = spot_size
            )
        }
    } else {
        if (is.null(colors)) {
            p <- vis_gene(
                spe_small,
                sampleid = sample_id, image_id = image_id, geneid = var_name,
                multi_gene_method = multi_gene_method, return_plots = TRUE,
                spatial = spatial, point_size = spot_size, assayname = assayname,
                alpha = NA, auto_crop = FALSE,
                minCount = minCount, na_color = NA
            )
        } else {
            p <- vis_gene(
                spe_small,
                sampleid = sample_id, image_id = image_id, geneid = var_name,
                multi_gene_method = multi_gene_method, return_plots = TRUE,
                spatial = spatial, point_size = spot_size, assayname = assayname,
                cont_colors = colors, alpha = NA, auto_crop = FALSE,
                minCount = minCount, na_color = NA
            )
        }
    }

    #   Remove the legend if requested
    if (!include_legend) {
        p <- p + theme(legend.position = "none")
    }

    #   Overwrite the title
    p <- p + labs(title = title)

    return(p)
}

plotDir <- here("plots", "14_pseudobulk_spatial", "01_precast", "plot_select_markers")
subregion_genes <- list(
    "shell" = c(
        "CALB2", "CPNE4", "ARHGAP6", "OPRK1", "HTR4"
    ),
    "core" = toupper(
        c(
            "CALB1", "LAMP5", "PEG10", "PDE1A", "ADORA2A"
        )
    ),
    "D1_islands" = c("FOXP2", "OPRM1", "RXFP1", "CPNE4", "DRD1", "PDYN", "PENK"), 
    "white_matter" = c("MBP", "MOBP", "PPP1R1B")
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
k <- 10
precast_path = here(
    'processed-data', '10_precast', 'nnSVG_precast', sprintf('PRECAST_k%s.csv', k)
)
cluster_col <- paste0("precast_k", k)
spe[[cluster_col]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(precast_path), by = 'key') |>
    pull(cluster) |>
    as.factor()

spe <- spe[ ,!is.na(spe[[cluster_col]])]
rownames(spe) <- rowData(spe)$gene_name

# Make preliminary cluster plots
pCluster <- list()
for(donor in select_samples){
    pCluster[[donor]] <- spot_plot(spe, sample_id = donor, var_name = paste0("precast_k", k),
 is_discrete = TRUE, spatial = TRUE, title = donor) + guides(fill = guide_legend(override.aes = list(size = 5)))  + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), title =element_text(size=12, face='bold'), 
 panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
}

pdf(file.path(plotDir, "cluster_assignments.pdf"))
pCluster[["Br6432"]]
pCluster[["Br6522"]]
pCluster[["Br3942"]]
dev.off()

for(i in c(1:length(subregion_genes$D1_islands))){
    cat(subregion_genes$D1_islands[i], "\n")
    p_D1islands <- list()
    for(donor in select_samples){
        p_D1islands[[donor]] <- spot_plot2(spe, sample_id = donor, var_name = subregion_genes$D1_islands[i], 
    is_discrete = FALSE, spatial = TRUE, minCount = 0.5) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), title =element_text(size=12, face='bold'), 
 panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
    }
    pdf(file.path(plotDir, "D1_islands", paste0(subregion_genes$D1_islands[i], ".pdf")), width = 6, height = 6)
    p_D1islands[["Br6432"]]
    p_D1islands[["Br6522"]]
    p_D1islands[["Br3942"]]
    dev.off()
}

for(i in c(1:length(subregion_genes$core))){
    cat(subregion_genes$core[i], "\n")
    p_core <- list()
    for(donor in select_samples){
        p_core[[donor]] <- spot_plot2(spe, sample_id = donor, var_name = subregion_genes$core[i], 
    is_discrete = FALSE, spatial = TRUE, colors = viridisLite::viridis(21), minCount = 0.5) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), title =element_text(size=12, face='bold'), 
 panel.border = element_rect(colour = "black", fill=NA, size=0.5))
    }
    pdf(file.path(plotDir, "core", paste0(subregion_genes$core[i], ".pdf")), width = 6, height = 6)
    p_core[["Br6432"]]
    p_core[["Br6522"]]
    p_core[["Br3942"]]
    dev.off()
}

for(i in c(1:length(subregion_genes$shell))){
    cat(subregion_genes$shell[i], "\n")
    p_shell <- list()
    for(donor in select_samples){
        p_shell[[donor]] <- spot_plot2(spe, sample_id = donor, var_name = subregion_genes$shell[i], 
    is_discrete = FALSE, spatial = TRUE, colors = viridisLite::viridis(21), minCount = 0.5) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), title =element_text(size=12, face='bold'), 
 panel.border = element_rect(colour = "black", fill=NA, size=0.5))
    }
    pdf(file.path(plotDir, "shell", paste0(subregion_genes$shell[i], ".pdf")), width = 6, height = 6)
    p_shell[["Br6432"]]
    p_shell[["Br6522"]]
    p_shell[["Br3942"]]
    dev.off()
}

for(i in c(1:length(subregion_genes$white_matter))){
    p_white_matter <- list()
    for(donor in select_samples){
        p_white_matter[[donor]] <- spot_plot2(spe, sample_id = donor, var_name = subregion_genes$white_matter[i], 
    is_discrete = FALSE, spatial = TRUE, colors = viridisLite::viridis(21), minCount = 0.5) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), title =element_text(size=12, face='bold'), 
 panel.border = element_rect(colour = "black", fill=NA, size=0.5))
    }
    pdf(file.path(plotDir, "white_matter", paste0(subregion_genes$white_matter[i], ".pdf")), width = 6, height = 6)
    p_white_matter[["Br6432"]]
    p_white_matter[["Br6522"]]
    p_white_matter[["Br3942"]]
    dev.off()
}