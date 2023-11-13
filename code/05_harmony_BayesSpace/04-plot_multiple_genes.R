library(SpatialExperiment)
library(here)
library(sessioninfo)
library(spatialLIBD)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(spatialNAcUtils)
library(HDF5Array)

################################################################################
#   Path and variable definitions
################################################################################

dlpfc_repo = '/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC'
spe_nac_dir = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered_hdf5'
)
spe_dlpfc_path = file.path(
    dlpfc_repo, 'processed-data/spot_deconvo/05-shared_utilities/IF/spe.rds'
)
marker_dlpfc_path = file.path(
    dlpfc_repo,
    'processed-data/spot_deconvo/05-shared_utilities/marker_stats_layer.rds'
)
plot_dir = here('plots', '05_harmony_BayesSpace', 'multi_gene')
z_score_dir = file.path(plot_dir, 'z_score')
sparsity_dir = file.path(plot_dir, 'sparsity')
pc_dir = file.path(plot_dir, 'pc')

num_markers = 25

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
    'white_matter' = c('MBP', 'GFAP', 'PLP1', 'AQP4'),
    'D1_islands' = c('OPRM1', 'CHST9')
)

best_sample_dlpfc = 'Br6522_Ant_IF'
best_sample_nac = 'Br8492'

plot_funs = c(spot_plot_z_score, spot_plot_sparsity, spot_plot_pca)
names(plot_funs) = c(z_score_dir, sparsity_dir, pc_dir)

################################################################################
#   Functions
################################################################################

#   spe: SpatialExperiment with rowData column 'gene_id' corresponding the
#       Ensembl gene names, and colData column 'sample_id'
#   genes: named list of character vectors whose names are "subregion"s for
#       which the genes (as Ensembl IDs) are markers
#   dataset: character(1) included in output PDF names
#   best_sample: character(1) in spe$sample_id giving a specific sample that
#       be plotted across each subregion for comparison
#   plot_dir: character(1) giving directory to write plots to
#   plot_fun: function from spatialNAcUtils for plotting multiple genes
#       spatially
#
#   For each subregion and sample, plot a spatial map summarizing
#   gene expression across 'genes', using the method described in 'plot_fun'
make_plots = function(
        spe, genes, dataset, best_sample, plot_dir, plot_fun
    ) {
    plot_list_sample = list()

    #   Plot all samples together in one PDF page, but have separate PDFs for
    #   each subregion
    for (subregion in sort(names(genes))) {
        plot_list_subregion = list()
        for (sample_id in unique(spe$sample_id)) {
            if (dataset == 'DLPFC') {
                plot_list_subregion[[sample_id]] = plot_fun(
                    spe, genes[[subregion]], sample_id, title = sample_id, 
                    include_legend = TRUE, assayname = "counts"
                )
            } else {
                plot_list_subregion[[sample_id]] = plot_fun(
                    spe, genes[[subregion]], sample_id, title = sample_id, 
                    include_legend = TRUE
                )
            }
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

spe_nac = loadHDF5SummarizedExperiment(spe_nac_dir)

#   Convert NAc gene symbols to Ensembl IDs. All gene symbols should exist in
#   the SpatialExperiment
stopifnot(all(unlist(genes_nac) %in% rowData(spe_nac)$gene_name))
for (subregion in names(genes_nac)) {
    genes_nac[[subregion]] = rowData(spe_nac)[
        match(genes_nac[[subregion]], rowData(spe_nac)$gene_name),
        'gene_id'
    ]
}

dir.create(z_score_dir, showWarnings = FALSE)
dir.create(sparsity_dir, showWarnings = FALSE)
dir.create(pc_dir, showWarnings = FALSE)

for (i in 1:length(plot_funs)) {
    make_plots(
        spe = spe_nac,
        genes = genes_nac,
        dataset = 'NAc',
        best_sample = best_sample_nac,
        plot_dir = names(plot_funs)[[i]],
        plot_fun = plot_funs[[i]]
    )
}

#-------------------------------------------------------------------------------
#   DLPFC plots
#-------------------------------------------------------------------------------

spe_dlpfc = readRDS(spe_dlpfc_path)
spe_dlpfc$exclude_overlapping = FALSE

#   Grab just the top [num_markers] markers for each excitatory cell type
#   marking for 1 specific layer
marker_dlpfc = readRDS(marker_dlpfc_path) |>
    mutate(
        symbol = rowData(spe_dlpfc)$gene_name[
            match(gene, rownames(spe_dlpfc))
        ]
    ) |>
    filter(
        #   Just one-layer excitatory cell types
        grepl('^Excit_L[3-6]$', cellType.target),
        #   No mitochondrial genes
        !grepl("^MT-", symbol)
    )

#   "Re-rank" rank_ratio, since there may be missing ranks now
for (ct in unique(marker_dlpfc$cellType.target)) {
    old_ranks <- marker_dlpfc |>
        filter(cellType.target == ct) |>
        pull(rank_ratio) |>
        sort()

    for (i in 1:length(which((marker_dlpfc$cellType.target == ct)))) {
        index <- which(
            (marker_dlpfc$cellType.target == ct) &
                (marker_dlpfc$rank_ratio == old_ranks[i])
        )
        stopifnot(length(index) == 1)
        marker_dlpfc[index, "rank_ratio"] <- i
    }
}

#   Now take the top [num_markers] markers
marker_dlpfc = marker_dlpfc |>
    filter(rank_ratio <= num_markers)

#   Form a list of genes where target cell types are the names and Ensembl IDs
#   are the values
genes_dlpfc = list()
for (cell_type in unique(marker_dlpfc$cellType.target)) {
    genes_dlpfc[[cell_type]] = marker_dlpfc |>
        filter(cellType.target == cell_type) |>
        pull(gene)
}

#   All genes (originally from 'marker_dlpfc') should be in 'spe_dlpfc'
stopifnot(all(unlist(genes_dlpfc) %in% rowData(spe_dlpfc)$gene_id))

for (i in 1:length(plot_funs)) {
    make_plots(
        spe = spe_dlpfc,
        genes = genes_dlpfc,
        dataset = 'DLPFC',
        best_sample = best_sample_dlpfc,
        plot_dir = names(plot_funs)[[i]],
        plot_fun = plot_funs[[i]]
    )
}

session_info()
