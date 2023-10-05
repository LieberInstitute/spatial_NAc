library(SpatialExperiment)
library(here)
library(sessioninfo)
library(spatialLIBD)
library(ggplot2)
library(cowplot)

source(here('code', '07_spot_deconvo', 'shared_functions.R'))

spe_path = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered.rds'
)
plot_dir = here('plots', '05_harmony_BayesSpace', 'multi_gene')
plot_nrow = 2

dir.create(plot_dir, showWarnings = FALSE)

spe = readRDS(spe_path)
spe$exclude_overlapping = spe$exclude_overlapping == "True"

#   List of genes provided by Svitlana, named by the subregion they're markers
#   for
genes = list(
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

#   All gene symbols should exist in the SpatialExperiment
stopifnot(all(unlist(genes) %in% rowData(spe)$gene_name))

#   Plot all samples together in one PDF page, but have separate PDFs for each
#   subregion
for (subregion in names(genes)) {
    #   Subset lognorm counts to just the selected genes
    gene_counts = assays(spe)$logcounts[
        match(genes[[subregion]], rowData(spe)$gene_name),
    ]

    #   For each spot, average expression Z-scores across all selected genes
    gene_z = (gene_counts - rowMeans(gene_counts)) / (rowSdDiffs(gene_counts))
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
    pdf(
        file.path(plot_dir, paste0(subregion, '.pdf')),
        width = 7 * ceiling(length(unique(spe$sample_id)) / plot_nrow),
        height = plot_nrow * 7
    )
    print(plot_grid(plotlist = plot_list, nrow = plot_nrow))
    dev.off()
}

session_inf