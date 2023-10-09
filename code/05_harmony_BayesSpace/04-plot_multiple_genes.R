library(SpatialExperiment)
library(here)
library(sessioninfo)
library(spatialLIBD)
library(ggplot2)
library(cowplot)

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

################################################################################
#   Functions
################################################################################

#   Given a SpatialExperiment and set of genes (contained in the 'gene_name'
#   rowData column), export a set of spatial plots that combine all genes'
#   expression (by subregion, expected to correspond to the names of 'genes').
#   The average (across genes) is taken of each gene's Z-score (across all spots
#   in all samples).
plot_z_score = function(spe, genes) {
    #   Plot all samples together in one PDF page, but have separate PDFs for
    #   each subregion
    for (subregion in names(genes)) {
        #   Subset lognorm counts to just the selected genes
        gene_counts = assays(spe)$logcounts[
            match(genes[[subregion]], rowData(spe)$gene_name),
        ]

        #   For each spot, average expression Z-scores across all selected genes
        gene_z = (gene_counts - rowMeans(gene_counts)) /
            (rowSdDiffs(gene_counts))
        spe$temp_var = colMeans(gene_z)

        #   Plot spatial distribution of averaged expression Z-scores for each
        #   sample (donor)
        plot_list = list()
        for (sample_id in unique(spe$sample_id)) {
            plot_list[[sample_id]] = spot_plot(
                spe, sample_id, title = sample_id, var_name = 'temp_var',
                include_legend = TRUE, is_discrete = FALSE, minCount = 0
            )
        }

        #   Save plots
        pdf(file.path(plot_dir, paste0(subregion, '.pdf')))
        print(plot_list)
        dev.off()
    }
}

################################################################################
#   Main
################################################################################

dir.create(plot_dir, showWarnings = FALSE)

spe_nac = readRDS(spe_path)
spe_nac$exclude_overlapping = spe_nac$exclude_overlapping == "True"

#   All gene symbols should exist in the SpatialExperiment for NAc
stopifnot(all(unlist(genes_nac) %in% rowData(spe_nac)$gene_name))

#   Plot NAc data using Z-score approach
plot_z_score(spe_nac, genes_nac)

session_info()
