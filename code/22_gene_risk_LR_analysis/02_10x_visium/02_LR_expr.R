library(spacexr)
library(Matrix)
library(SingleCellExperiment)
library(here)
library(scran)
library(scater)
library(SpatialExperiment)
library(spatialLIBD)
library(spatialNAcUtils)
library(HDF5Array)
library(ggplot2)
library(patchwork)
library(getopt)
library(viridis)
library(ggpubr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(glue)
codeDir <- here("code")
source(file.path(codeDir, "plot_utils.R"))

# Read in SPE, add cell type deconvolution, and PRECAST clusters
opt <- list()
opt$marker_genes <- TRUE

spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5")
spe <- loadHDF5SummarizedExperiment(spe_dir)
spe_Br2743 <- mirrorObject(spe, sample_id = "Br2743", image_id = "lowres", axis = "v")
spe_Br8492 <- mirrorObject(spe, sample_id = "Br8492", image_id = "lowres", axis = "v")
spe_Br8325 <- mirrorObject(spe, sample_id = "Br8325", image_id = "lowres", axis = "v")
spe_Br3942 <- mirrorObject(spe, sample_id = "Br3942", image_id = "lowres", axis = "v")

spe <- spe[ ,!spe$sample_id %in% c("Br2743", "Br8492", "Br8325", "Br3942")]
spe <- cbind(spe, spe_Br2743)
spe <- cbind(spe, spe_Br8492)
spe <- cbind(spe, spe_Br8325)
spe <- cbind(spe, spe_Br3942)

sample_ids <- levels(spe$sample_id)
RCTD_list <- lapply(sample_ids, function(isample){
    RCTD_dir <- here::here("processed-data", "08_spot_deconvo", "01_RCTD", isample)
    if(opt$marker_genes){
        myRCTD <- readRDS(file.path(RCTD_dir, "results_RCTD_markers.rds"))
    }else{
        myRCTD <- readRDS(file.path(RCTD_dir, "results_RCTD.rds"))
    }
    myRCTD
})

weights <- lapply(RCTD_list, function(myRCTD){
    my_results = myRCTD@results
    my_weights = lapply(my_results, function(x) x$all_weights)
    my_weights_df <- data.frame(do.call(rbind, my_weights))
    my_weights_df
})
weights <- data.frame(do.call(rbind, weights))

coords <- lapply(RCTD_list, function(myRCTD){
    my_coords = myRCTD@spatialRNA@coords
    my_coords
})
coords <- data.frame(do.call(rbind, coords))

spe <- spe[ ,colnames(spe) %in% rownames(coords)]
spe <- spe[ ,match(rownames(coords), colnames(spe))]
colData(spe) <- cbind(colData(spe), weights)

# Add the final clusters
# Add precast results to the spe object
clusters_file <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters", "precast_clusters.csv")
spe[["precast_clusters"]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(clusters_file), by = 'key') |>
    pull(cluster) |>
    as.factor()
# Remove spots with no PRECAST output
spe <- spe[ ,!is.na(spe[["precast_clusters"]])]

# Specify the colors for the domains
safe_colorblind_palette <- c("#66A61E","#1B9E77", "#7570B3","#E7298A","#D95F02" , "#E6AB02","#666666", "#A6761D")
spe$precast_clusters <- factor(spe$precast_clusters, 
levels = c("MSN 1", "MSN 2", "MSN 3", "D1 islands", "Excitatory", "Inhibitory", "WM", "Endothelial/Ependymal"))

rownames(spe) <- make.unique(rowData(spe)$gene_name)

plot_genes_expr <- function(spe, genes, cluster_col = "precast_clusters"){
    expr_mat <- assays(spe)$logcounts
    expr_sub <- as.matrix(t(expr_mat[genes, ]))
    plot_df <- data.frame(expr_sub)
    plot_df$Domains <- spe[[cluster_col]]

    plot_df <- reshape2::melt(plot_df, id.vars = "Domains")
    colnames(plot_df) <- c("Domains", "Gene", "Expr")
    levels(plot_df$Gene) <- genes

    cat_colors <- c("#d73027", "#91bfdb", "#4575b4")
    p <- ggplot(plot_df, aes(x = Domains, y = Expr)) +
    geom_boxplot(aes(fill = Gene), color = "black", outliers = FALSE)+
    facet_wrap(~ Gene, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = cat_colors) + 
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey30"),
      strip.text = element_text(color = "white", face = "bold")
    ) +
    labs(
      title = "",
      x = NULL,
      y = "Expression"
    )

  return(p)
}


p1 <- plot_genes_expr(spe, c("OPRM1", "PENK", "PDYN"))

plot_dir <- here("plots", "22_gene_risk_LR_analysis", "05-10x_visium")
pdf(file.path(plot_dir, "LR_pairs_expression_OPRM1_PENK_PDYN.pdf"), width = 5, height = 6)
print(p1)
dev.off()


plot_cell_types <- function(spe, cell_types, cluster_col = "precast_clusters"){
    deconv <- colData(spe)[ ,cell_types]
    plot_df <- data.frame(deconv)
    plot_df$Domains <- spe[[cluster_col]]
    plot_df <- reshape2::melt(plot_df, id.vars = "Domains")
    colnames(plot_df) <- c("Domains", "Cell_type", "value")
    levels(plot_df$Cell_type) <- cell_types

    cell_type_colors <-  c(
    "Astrocyte_A" = "#D95F00",
    "Astrocyte_B" = "#FF0016",
    "DRD1_MSN_A" = "#ECA31C",
    "DRD1_MSN_B" = "#D079AA",
    "DRD1_MSN_C" = "#352391",
    "DRD1_MSN_D" = "#d156aa",
    "DRD2_MSN_A" = "#58B6ED",
    "DRD2_MSN_B" ="#F80091",
    "Endothelial" = "#D00DFF",
    "Ependymal" = "#0D73B4",
    "Excitatory" = "#5C6300",
    "Inh_A" = "#35FB00",
    "Inh_B" = "#7A0096",
    "Inh_C" = "#854222", 
    "Inh_D" = "#A7F281", 
    "Inh_E" = "#0DFBFA", 
    "Inh_F" = "black",
    "Microglia" = "#F2E642",
    "Oligo" = "#4F4753",
    "OPC" = "#0D9F72")

    p <- ggplot(plot_df, aes(x = Domains, y = value)) +
    geom_boxplot(aes(fill = Cell_type), color = "black", outliers = FALSE)+
    facet_wrap(~ Cell_type, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = cell_type_colors) + 
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey30"),
      strip.text = element_text(color = "white", face = "bold")
    ) +
    labs(
      title = "",
      x = NULL,
      y = "RCTD scores"
    )

  return(p)
}


p2 <- plot_cell_types(spe, c("DRD2_MSN_A", "DRD1_MSN_A", "DRD1_MSN_B","DRD1_MSN_D"))

pdf(file.path(plot_dir, "LR_pairs_cell_types_OPRM1.pdf"), width = 7, height = 5)
print(p2)
dev.off()

make_gene_escheR <- function(spe, genes, cluster_col = "precast_clusters"){
  sample_order <- c("Br2743", "Br6432", "Br6423", "Br2720", "Br6471", "Br6522","Br8492", "Br8325", "Br8667", "Br3942")
  spe$precast_clusters <- factor(as.character(spe$precast_clusters), levels = c("D1 islands", "Endothelial/Ependymal", "Excitatory", "Inhibitory", "MSN 1", "MSN 2", "MSN 3", "WM")) 
  safe_colorblind_palette <- c("#E7298A", "#A6761D","#D95F02" , "#E6AB02",  "#66A61E","#1B9E77", "#7570B3","#666666")
  names(safe_colorblind_palette) <- levels(spe$precast_clusters)
  for (gene in genes) {
    print(gene)
    plot_list <- list()
    for (donor in sample_order) {
      spe_sub <- spe[ ,spe$sample_id == donor]
      spe_sub <- spe_sub[ ,!spe_sub$exclude_overlapping]
      expr_mat <- assays(spe_sub)$logcounts
      spe_sub$gene_expr <- as.numeric(expr_mat[gene, ])
      plot_list[[donor]] <- make_escheR(spe_sub) |> 
                          add_fill(var = "gene_expr") |> 
                          add_ground(var = "precast_clusters", stroke = 0.5) +
                          scale_color_manual(values = safe_colorblind_palette) +
                          scale_fill_gradient(low = "white", high = "black")
    }
  pdf(file.path(plot_dir, paste0(gene, "_escheR.pdf")), width = 12, height = 12)
  for (p in plot_list) print(p)
  dev.off()
 }
 }

make_gene_escheR(spe, c("OPRM1", "PENK", "PDYN"))