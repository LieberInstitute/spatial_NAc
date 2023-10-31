library(SpatialExperiment)
library(here)
library(sessioninfo)
library(spatialLIBD)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(spatialNAcUtils)

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
    'white_matter' = c('MBP', 'GFAP', 'PLP1', 'AQP4')
)

best_sample_dlpfc = 'Br6522_Ant_IF'
best_sample_nac = 'Br8492'

################################################################################
#   Functions
################################################################################

#   spe: SpatialExperiment with rowData column 'gene_id' corresponding the
#       Ensembl gene names, and colData column 'sample_id'
#   genes: named list of character vectors whose names are "subregion"s for
#       which the genes (as Ensembl IDs) are markers
#   dataset: character(1) included in output PDF names
#   assay: character(1) expected to be in names(assays(spe))
#   best_sample: character(1) in spe$sample_id giving a specific sample that
#       be plotted across each subregion for comparison
#   plot_dir: character(1) giving directory to write plots to
#
#   For each subregion and sample, plot a spatial map of the average
#   (across genes) of each gene's Z-score of expression (across all spots
#   in all samples). Returns NULL
plot_z_score = function(
        spe, genes, dataset, assay, best_sample, plot_dir
    ) {
    plot_list_sample = list()

    #   Plot all samples together in one PDF page, but have separate PDFs for
    #   each subregion
    for (subregion in sort(names(genes))) {
        #   Subset assay to just the selected genes
        gene_counts = assays(spe)[[assay]][
            match(genes[[subregion]], rowData(spe)$gene_id),
        ]

        #   For each spot, average expression Z-scores across all selected genes
        gene_z = (gene_counts - rowMeans(gene_counts)) /
            (rowSdDiffs(gene_counts))
        spe$temp_var = colMeans(gene_z, na.rm = TRUE)

        #   Plot spatial distribution of averaged expression Z-scores for each
        #   sample
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

#   spe: SpatialExperiment with counts assay and rowData column 'gene_id'
#       corresponding the Ensembl gene names, and colData column 'sample_id'
#   genes: named list of character vectors whose names are "subregion"s for
#       which the genes (as Ensembl IDs) are markers
#   dataset: character(1) included in output PDF names
#   best_sample: character(1) in spe$sample_id giving a specific sample that
#       be plotted across each subregion for comparison
#   plot_dir: character(1) giving directory to write plots to
#
#   For each subregion and sample, plot a spatial map of the proportion of
#   markers for the given subregion that have nonzero expression (computed in
#   each spot). Returns NULL
plot_sparsity = function(spe, genes, dataset, best_sample, plot_dir) {
    plot_list_sample = list()

    #   Plot all samples together in one PDF page, but have separate PDFs for
    #   each subregion
    for (subregion in sort(names(genes))) {
        #   For each spot, compute proportion of marker genes with nonzero
        #   expression
        spe$prop_nonzero_marker <- colMeans(
            assays(
                    spe[match(genes[[subregion]], rowData(spe)$gene_id),]
                )$counts > 0
        )

        #   Plot spatial distribution of this proportion for each sample (donor)
        plot_list_subregion = list()
        for (sample_id in unique(spe$sample_id)) {
            plot_list_subregion[[sample_id]] = spot_plot(
                spe, sample_id, title = sample_id,
                var_name = 'prop_nonzero_marker', include_legend = TRUE,
                is_discrete = FALSE, minCount = 0.1
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

#   spe: SpatialExperiment with rowData column 'gene_id' corresponding the
#       Ensembl gene names, and colData column 'sample_id'
#   genes: named list of character vectors whose names are "subregion"s for
#       which the genes (as Ensembl IDs) are markers
#   dataset: character(1) included in output PDF names
#   assay: character(1) expected to be in names(assays(spe))
#   best_sample: character(1) in spe$sample_id giving a specific sample that
#       be plotted across each subregion for comparison
#   plot_dir: character(1) giving directory to write plots to
#
#   For each subregion and sample, plot a spatial map summarizing expression of
#   all 'genes'. Gene expression is subsetted to the genes of interest, and PCA
#   is performed on the transposed expression data, such that each PC represents
#   a "spot profile", where larger elements represent spots of greater
#   expression variation within the gene set. Returns NULL
plot_pc = function(spe, genes, dataset, assay, best_sample, plot_dir) {
    plot_list_sample = list()

    #   Plot all samples together in one PDF page, but have separate PDFs for
    #   each subregion
    for (subregion in sort(names(genes))) {
        #   Take gene expression for the select genes, and transpose it. Then
        #   perform PCA, so PCs capture highly variable spots with respect to
        #   the marker genes
        spot_exp = t(assays(spe)[[assay]][genes[[subregion]],])
        pc_exp = prcomp(spot_exp, center = TRUE, scale = TRUE)
        spe$pc_select_genes <- pc_exp$x[,'PC1']

        #   Given that:
        #       - 'genes' is assumed to represent markers of the subregion (and
        #         thus their expression follows a similar pattern spatially)
        #       - the first PC captures the most variation
        #   Then each gene's coefficients to the first PC should tend to have
        #   the same sign. Next, the sign of each PC is arbitary, and we'd like
        #   plots to have positive values where expression is greater. If most
        #   genes have negative coefficients to the first PC, we reverse the
        #   sign of the coefficients to make visual intrepretation consistent
        if (mean(pc_exp$rotation[,1] > 0) < 0.5) {
            spe$pc_select_genes = -1 * spe$pc_select_genes
        }

        #   Plot spatial distribution of this proportion for each sample (donor)
        plot_list_subregion = list()
        for (sample_id in unique(spe$sample_id)) {
            plot_list_subregion[[sample_id]] = spot_plot(
                spe, sample_id, title = sample_id,
                var_name = 'pc_select_genes', include_legend = TRUE,
                is_discrete = FALSE, minCount = 0
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

z_score_dir = file.path(plot_dir, 'z_score')
sparsity_dir = file.path(plot_dir, 'sparsity')
pc_dir = file.path(plot_dir, 'pc')
dir.create(z_score_dir, showWarnings = FALSE)
dir.create(sparsity_dir, showWarnings = FALSE)
dir.create(pc_dir, showWarnings = FALSE)

#   Plot NAc data using Z-score approach
plot_z_score(
    spe = spe_nac,
    genes = genes_nac,
    dataset = 'NAc',
    assay = 'logcounts',
    best_sample = best_sample_nac,
    plot_dir = z_score_dir
)

#   Plot NAc data using sparsity approach
plot_sparsity(
    spe = spe_nac,
    genes = genes_nac,
    dataset = 'NAc',
    best_sample = best_sample_nac,
    plot_dir = sparsity_dir
)

plot_pc(
    spe = spe_nac,
    genes = genes_nac,
    dataset = 'NAc',
    assay = 'logcounts',
    best_sample = best_sample_nac,
    plot_dir = pc_dir
)

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

#   Plot DLPFC data using Z-score approach
plot_z_score(
    spe = spe_dlpfc,
    genes = genes_dlpfc,
    dataset = 'DLPFC',
    assay = 'counts',
    best_sample = best_sample_dlpfc,
    plot_dir = z_score_dir
)

#   Plot DLPFC data using sparsity approach
plot_sparsity(
    spe = spe_dlpfc,
    genes = genes_dlpfc,
    dataset = 'DLPFC',
    best_sample = best_sample_dlpfc,
    plot_dir = sparsity_dir
)

session_info()
