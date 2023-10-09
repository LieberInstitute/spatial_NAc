library(SpatialExperiment)
library(here)
library(sessioninfo)
library(spatialLIBD)
library(ggplot2)
library(cowplot)
library(tidyverse)

source(here('code', '07_spot_deconvo', 'shared_functions.R'))

################################################################################
#   Path and variable definitions
################################################################################

dlpfc_repo = '/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC'
spe_nac_path = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered.rds'
)
spe_dlpfc_path = file.path(
    dlpfc_repo, 'processed-data/spot_deconvo/05-shared_utilities/IF/spe.rds'
)
marker_dlpfc_path = file.path(
    dlpfc_repo,
    'processed-data/spot_deconvo/05-shared_utilities/marker_stats_layer.rds'
)
plot_dir = here('plots', '05_harmony_BayesSpace', 'multi_gene')

num_markers = 10

#   List of genes provided by Svitlana, named by the subregion they're markers
#   for
genes_nac = list(
    'shell' = c(
        'TAC3', 'TAC1', 'CALB2', 'GRIA4', 'CPNE4', 'GREB1L', 'ARHGAP6', 'OPRK1'
    ),
    'core' = toupper(
        c(
            'CALB1', 'HPCAL4', 'ZDBF2', 'LAMP5', 'Scn4b', 'Rgs9', 'Ppp3ca',
            'Oprd1', 'Adora2a', 'Pde1a', 'Peg10', 'Dlk1', 'Htr2c'
        )
    ),
    'white_matter' = c('MBP', 'GFAP', 'PLP1', 'AQP4')
)

best_sample_dlpfc = 'Br6522_Ant_IF'
best_sample_nac = 'Br8492'

################################################################################
#   Functions
################################################################################

#   Given a SpatialExperiment and set of genes (contained in the 'gene_id'
#   rowData column), export a set of spatial plots that combine all genes'
#   expression (by subregion, expected to correspond to the names of 'genes').
#   The average (across genes) is taken of each gene's Z-score (across all spots
#   in all samples).
plot_z_score = function(spe, genes, dataset, assay, best_sample) {
    plot_list_sample = list()

    #   Plot all samples together in one PDF page, but have separate PDFs for
    #   each subregion
    for (subregion in sort(names(genes))) {
        #   Subset lognorm counts to just the selected genes
        gene_counts = assays(spe)[[assay]][
            match(genes[[subregion]], rowData(spe)$gene_id),
        ]

        #   For each spot, average expression Z-scores across all selected genes
        gene_z = (gene_counts - rowMeans(gene_counts)) /
            (rowSdDiffs(gene_counts))
        spe$temp_var = colMeans(gene_z)
        if(any(is.na(spe$temp_var))) {
            stop('some unexpressed genes')
        }

        #   Plot spatial distribution of averaged expression Z-scores for each
        #   sample (donor)
        plot_list_subregion = list()
        for (sample_id in unique(spe$sample_id)) {
            plot_list_subregion[[sample_id]] = spot_plot(
                spe, sample_id, title = sample_id, var_name = 'temp_var',
                include_legend = TRUE, is_discrete = FALSE, minCount = 0,
                assayname = assay
            )
        }

        plot_list_sample[[subregion]] = plot_list_subregion[[best_sample]] +
            labs(title = subregion)

        #   Save plots
        pdf(file.path(plot_dir, sprintf('%s_%s.pdf', dataset, subregion)))
        print(plot_list_subregion)
        dev.off()
    }

    #   Now plot every subregion for the best sample (a single-page, single-row
    #   PDF)
    pdf(
        file.path(plot_dir, sprintf('%s_%s.pdf', dataset, best_sample)),
        width = 7 * length(genes)
    )
    print(plot_grid(plotlist = plot_list_sample, nrow = 1))
    dev.off()
}

################################################################################
#   Main
################################################################################

dir.create(plot_dir, showWarnings = FALSE)

#-------------------------------------------------------------------------------
#   NAc plots
#-------------------------------------------------------------------------------

spe_nac = readRDS(spe_nac_path)
spe_nac$exclude_overlapping = spe_nac$exclude_overlapping == "True"

#   Convert NAc gene symbols to Ensembl IDs. All gene symbols should exist in
#   the SpatialExperiment
stopifnot(all(unlist(genes_nac) %in% rowData(spe_nac)$gene_name))
for (subregion in names(genes_nac)) {
    genes_nac[[subregion]] = rowData(spe_nac)[
        match(genes_nac[[subregion]], rowData(spe_nac)$gene_name),
        'gene_id'
    ]
}

#   Plot NAc data using Z-score approach
plot_z_score(spe_nac, genes_nac, 'NAc', best_sample_nac)

#-------------------------------------------------------------------------------
#   DLPFC plots
#-------------------------------------------------------------------------------

spe_dlpfc = readRDS(spe_dlpfc_path)
spe_dlpfc$exclude_overlapping = FALSE

#   Grab just the top [num_markers] markers for each excitatory cell type
#   marking for 1 specific layer
marker_dlpfc = readRDS(marker_dlpfc_path) |>
    filter(
        rank_ratio <= num_markers,
        grepl('^Excit_L[3-6]$', cellType.target)
    )

#   Form a list of genes where target cell types are the names and Ensembl IDs
#   are the values
genes_dlpfc = list()
for (cell_type in unique(marker_dlpfc$cellType.target)) {
    genes_dlpfc[[cell_type]] = marker_dlpfc |>
        filter(cellType.target == cell_type) |>
        pull(gene)
}

#   Plot DLPFC data using Z-score approach
plot_z_score(spe_dlpfc, genes_dlpfc, 'DLPFC', 'counts', best_sample_dlpfc)

session_info()
