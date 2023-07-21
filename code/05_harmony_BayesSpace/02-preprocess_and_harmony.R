## Required libraries
library("here")
library("SpatialExperiment")
library("spatialLIBD")
library("sessioninfo")
library("scran") ## requires uwot for UMAP
library("uwot")
library("scater")
library("BiocParallel")
library("PCAtools")
library("ggplot2")
library("Polychrome")
library("harmony")

tsne_perplex_vals = c("05", "20", "50", "80")

## Create output directories
dir_plots <- here("plots", "05_harmony_BayesSpace")
dir_rdata <- here("processed-data", "05_harmony_BayesSpace")
spe_in <- file.path(dir_rdata, "spe.rds")

dir.create(dir_plots, showWarnings = FALSE)
dir.create(dir_rdata, showWarnings = FALSE)
dir.create(file.path(dir_rdata, "clusters_graphbased"), showWarnings = FALSE)
dir.create(file.path(dir_rdata, "clusters_graphbased_cut_at"), showWarnings = FALSE)

set.seed(20230712)

## Load the data
spe <- readRDS(spe_in)

#   For testing
# spe = spe[sample(1:nrow(spe), 2000), sample(1:ncol(spe), 5000)]
# spe <- spe[
#     rowSums(assays(spe)$counts) > 0, colSums(assays(spe)$counts) > 0
# ]

message("Running quickCluster()")

Sys.time()
spe$scran_quick_cluster <- quickCluster(
    spe,
    BPPARAM = MulticoreParam(4),
    block = spe$sample_id,
    block.BPPARAM = MulticoreParam(4)
)
Sys.time()

message("Running computeSumFactors()")
Sys.time()
## Might be needed:
# options(error = recover)
spe <-
    computeSumFactors(spe,
        clusters = spe$scran_quick_cluster,
        BPPARAM = MulticoreParam(4)
    )
Sys.time()

table(spe$scran_quick_cluster)

message("Running checking sizeFactors()")
summary(sizeFactors(spe))

message("Running logNormCounts()")
spe <- logNormCounts(spe)

message("Running modelGeneVar()")
## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
dec <- modelGeneVar(spe,
    block = spe$sample_id,
    BPPARAM = MulticoreParam(4)
)

pdf(file.path(dir_plots, "scran_modelGeneVar.pdf"), useDingbats = FALSE)
mapply(function(block, blockname) {
    plot(
        block$mean,
        block$total,
        xlab = "Mean log-expression",
        ylab = "Variance",
        main = blockname
    )
    curve(metadata(block)$trend(x),
        col = "blue",
        add = TRUE
    )
}, dec$per.block, names(dec$per.block))
dev.off()

message("Running getTopHVGs()")
top.hvgs <- getTopHVGs(dec, prop = 0.1)
length(top.hvgs)

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)

save(top.hvgs,
    top.hvgs.fdr5,
    top.hvgs.fdr1,
    file = file.path(dir_rdata, "top.hvgs.Rdata")
)


message("Running runPCA()")
Sys.time()
spe <- runPCA(spe, subset_row = top.hvgs, ncomponents = 50)
Sys.time()

# make elbow plot to determine PCs to use
percent.var <- attr(reducedDim(spe, "PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow

pdf(
    file.path(dir_plots, "pca_elbow.pdf"),
    useDingbats = FALSE
)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()

## Perform harmony batch correction
message("Running RunHarmony()")
Sys.time()
spe <- RunHarmony(spe, "sample_id", verbose = FALSE)
Sys.time()

#   Perform dimensionality reduction using both PCA and harmony's reduced
#   dimensions
for (dimred_var in c("PCA", "HARMONY")) {
    #   Run TSNE with several perplexity values
    for (perplex in tsne_perplex_vals) {
        message(
            sprintf(
                "Running runTSNE() perplexity %s on %s dimensions",
                perplex, dimred_var
            )
        )
        Sys.time()
        spe <-
            runTSNE(
                spe,
                dimred = dimred_var,
                name = sprintf("TSNE_perplexity%s.%s", perplex, dimred_var),
                perplexity = as.integer(perplex)
            )
        Sys.time()
    }

    #   Explore TSNE results via plots
    pdf(
        file = file.path(
            dir_plots,
            sprintf("tSNE_perplexity%s_%s_sample_id.pdf", perplex, dimred_var)
        ),
        width = 9
    )
    ggplot(
            data.frame(
                reducedDim(
                    spe, sprintf("TSNE_perplexity%s.%s", perplex, dimred_var)
                )
            ),
            aes(x = TSNE1, y = TSNE2, color = factor(spe$sample_id))
        ) +
        geom_point() +
        labs(color = "sample_id") +
        theme_bw()
    dev.off()

    #   Also run UMAP
    message(sprintf("Running runUMAP() on %s dimensions", dimred_var))
    Sys.time()
    spe <- runUMAP(
        spe, dimred = dimred_var, name = sprintf("UMAP.%s", dimred_var),
        BPPARAM = MulticoreParam(4)
    )
    Sys.time()

    #   Explore UMAP results, coloring by both donor and sample ID
    for (color_var in c("sample_id", "donor")) {
        pdf(
            file = file.path(
                dir_plots, sprintf("UMAP_%s_%s.pdf", color_var, dimred_var)
            )
        )
        ggplot(
            data.frame(reducedDim(spe, sprintf("UMAP.%s", dimred_var))),
            aes(x = UMAP1, y = UMAP2, color = factor(spe[[color_var]]))
        ) +
            geom_point() +
            labs(color = color_var) +
            theme_bw()
        dev.off()
    }
}

## Perform graph-based clustering on batch corrected-data
message("Running buildSNNGraph() on HARMONY dimensions")
Sys.time()
g_k10 <- buildSNNGraph(spe, k = 10, use.dimred = "HARMONY")
Sys.time()
save(g_k10, file = file.path(dir_rdata, "g_k10_harmony.Rdata"))

message("Running cluster_walktrap()")
Sys.time()
g_walk_k10 <- igraph::cluster_walktrap(g_k10)
Sys.time()
save(g_walk_k10, file = file.path(dir_rdata, "g_walk_k10_harmony.Rdata"))

clust_k10 <- sort_clusters(g_walk_k10$membership)
spe$SNN_k10 <- clust_k10 ## Add this one to the SPE too
## Export for later use
cluster_export(
    spe,
    "SNN_k10",
    cluster_dir = file.path(dir_rdata, "clusters_graphbased"),
)

message("Running cut_at() from k = 4 to 28")
clust_k5_list <- lapply(4:28, function(n) {
    message(paste(Sys.time(), "n =", n))
    sort_clusters(igraph::cut_at(g_walk_k10, n = n))
})
names(clust_k5_list) <- paste0("SNN_k10_k", 4:28)

## Add clusters to spe colData
for (i in seq_along(names(clust_k5_list))) {
    colData(spe) <- cbind(colData(spe), clust_k5_list[i])
    ## Add proper name
    colnames(colData(spe))[ncol(colData(spe))] <- names(clust_k5_list)[i]

    ## Export for later use outside the SPE object
    cluster_export(
        spe,
        names(clust_k5_list)[i],
        cluster_dir = file.path(dir_rdata, "clusters_graphbased_cut_at")
    )
}

## make plot
sample_ids <- unique(colData(spe)$sample_id)
pdf(file = file.path(dir_plots, "graph_based_harmony.pdf"))
for (i in seq_along(sample_ids)) {
    for (j in seq_along(names(clust_k5_list))) {
        clus_vals <- unique(clust_k5_list[[j]])
        cols <- Polychrome::palette36.colors(length(clus_vals))
        names(cols) <- clus_vals

        my_plot <- vis_clus(
            spe = spe,
            clustervar = names(clust_k5_list)[j],
            sampleid = sample_ids[i],
            colors = cols,
            auto_crop = FALSE,
            assayname = 'counts',
            ... = paste0(" ", names(clust_k5_list)[j])
        )
        print(my_plot)
    }
}
dev.off()


#   Do offset so we can run BayesSpace. Not here that 'array_row' is not
#   constrained to have max value 77; we instead find the largest 'array_row'
#   value of any sample, and use it to ensure samples are at least 5 rows apart
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * (max(spe$array_row) + 5)
names(auto_offset_row) <- unique(spe$sample_id)
spe$row <- spe$array_row + auto_offset_row[spe$sample_id]
spe$col <- spe$array_col

## Save new SPE object
saveRDS(spe, file.path(dir_rdata, "spe_harmony.rds"))

## Object size in GB
## (do this near the end in case lobstr crashes, it's happened to me once)
lobstr::obj_size(spe)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
