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
library(escheR)
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
#### Make plots showing domain specificity of LR expression
plot_LR_expression_by_cluster_faceted <- function(
  spe,
  ligand,
  receptor,
  cluster_col = "precast_clusters",
  sample_col = "sample_id",
  min_expr = 0
) {
  expr_mat <- assays(spe)$logcounts

  # Check genes exist
  if (!(ligand %in% rownames(expr_mat)) | !(receptor %in% rownames(expr_mat))) {
    warning(paste("Missing gene:", ligand, "or", receptor))
    return(NULL)
  }

  # Binary expression
  lig_expr <- expr_mat[ligand, ] > min_expr
  rec_expr <- expr_mat[receptor, ] > min_expr

  # Data frame with cluster/sample/category
  df <- data.frame(
    sample_id = colData(spe)[[sample_col]],
    cluster = colData(spe)[[cluster_col]],
    coexpr = lig_expr & rec_expr,
    lig = lig_expr,
    rec = rec_expr
  ) |>
    filter(!is.na(cluster)) |>
    pivot_longer(
      cols = c("coexpr", "lig", "rec"),
      names_to = "Category",
      values_to = "Expressing"
    ) |>
    group_by(sample_id, cluster, Category) |>
    summarise(Proportion = mean(Expressing), .groups = "drop") |>
    mutate(Category = factor(Category,
                             levels = c("coexpr", "lig", "rec"),
                             labels = c(
                               paste0(ligand, " & ", receptor),
                               paste0(ligand),
                               paste0(receptor)
                             )))

  # Plot
  cat_labels <- c(
  paste0(ligand, " & ", receptor),
  paste0(ligand),
  paste0(receptor))

  cat_colors <- setNames(
  c("#d73027", "#91bfdb", "#4575b4"), cat_labels)
  p <- ggplot(df, aes(x = cluster, y = Proportion)) +
    geom_boxplot(aes(fill = Category), color = "black", outlier.size = 1)+
    facet_wrap(~ Category, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = cat_colors) + 
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey30"),
      strip.text = element_text(color = "white", face = "bold")
    ) +
    labs(
      title = paste(ligand, "&", receptor),
      x = NULL,
      y = "Proportion of Spots"
    )

  return(p)
}

lr_pairs <- list(
    c("SEMA3E", "PLXND1"),
    c("EFNA5", "EPHA5"),
    c("PENK", "OPRM1"), 
    c("PDYN", "OPRM1"),
    c("CRH", "CRHR1"), 
    c("TAC1", "TACR3"), 
    c("BDNF", "NTRK2"), 
    c("NPY", "MC4R"), 
    c("NPY", "DPP4"),
    c("NLGN3", "NRXN1"), 
    c("NLGN1", "NRXN1"),
    c("NLGN1", "NRXN3"), 
    c("CNTN3", "PTPRG"), 
    c("CNTN4", "PTPRG"), 
    c("CNTN6", "PTPRG")
)

plot_list <- lapply(lr_pairs, function(pair) {
    ligand <- pair[1]
    receptor <- pair[2]
    plot_LR_expression_by_cluster_faceted(spe, ligand, receptor)
})

plot_dir <- here("plots", "22_gene_risk_LR_analysis", "05-10x_visium")
pdf(file.path(plot_dir, "LR_pairs_coexpression_by_domain.pdf"))
for(i in c(1:length(plot_list))){
    print(plot_list[[i]])
}
dev.off()

###########################################################################
# 2. Make spot plots
###########################################################################
add_LR_expression_category <- function(spe, ligand, receptor, min_expr = 1, new_colname = NULL) {
  expr_mat <- assays(spe)$logcounts

  if (is.null(new_colname)) {
    new_colname <- paste0("LR_", ligand, "_", receptor)
  }

  if (!(ligand %in% rownames(expr_mat)) | !(receptor %in% rownames(expr_mat))) {
    warning(paste("Missing gene:", ligand, "or", receptor))
    colData(spe)[[new_colname]] <- NA
    return(spe)
  }

  lig_expr <- expr_mat[ligand, ] > min_expr
  rec_expr <- expr_mat[receptor, ] > min_expr

  expr_class <- case_when(
    lig_expr & rec_expr ~ "Coexpressing",
    lig_expr & !rec_expr ~ paste0(ligand, " only"),
    !lig_expr & rec_expr ~ paste0(receptor, " only"),
    !lig_expr & !rec_expr ~ "None"
  )

  colData(spe)[[new_colname]] <- factor(expr_class, levels = c("Coexpressing", paste0(ligand, " only"), paste0(receptor, " only"), "None"))

  return(spe)
}

spe <- add_LR_expression_category(spe, ligand = "EFNA5", receptor = "EPHA5")
spe <- add_LR_expression_category(spe, ligand = "SEMA3E", receptor = "PLXND1")
spe <- add_LR_expression_category(spe, ligand = "TAC1", receptor = "TACR3")
spe <- add_LR_expression_category(spe, ligand = "PENK", receptor = "OPRM1")
spe <- add_LR_expression_category(spe, ligand = "PDYN", receptor = "OPRM1")
spe <- add_LR_expression_category(spe, ligand = "NPY", receptor = "MC4R")
spe <- add_LR_expression_category(spe, ligand = "NPY", receptor = "DPP4")
spe <- add_LR_expression_category(spe, ligand = "BDNF", receptor = "NTRK2")
spe <- add_LR_expression_category(spe, ligand = "CRH", receptor = "CRHR1")
spe <- add_LR_expression_category(spe, ligand = "NLGN3", receptor = "NRXN1")
spe <- add_LR_expression_category(spe, ligand = "NLGN1", receptor = "NRXN1")
spe <- add_LR_expression_category(spe, ligand = "NLGN1", receptor = "NRXN3")
spe <- add_LR_expression_category(spe, ligand = "CNTN3", receptor = "PTPRG")
spe <- add_LR_expression_category(spe, ligand = "CNTN4", receptor = "PTPRG")
spe <- add_LR_expression_category(spe, ligand = "CNTN6", receptor = "PTPRG")


sample_order <- c("Br2743", "Br6432", "Br6423", "Br2720", "Br6471", "Br6522","Br8492", "Br8325", "Br8667", "Br3942")
spe$precast_clusters <- factor(as.character(spe$precast_clusters), levels = c("D1 islands", "Endothelial/Ependymal", "Excitatory", "Inhibitory", "MSN 1", "MSN 2", "MSN 3", "WM")) 
safe_colorblind_palette <- c("#E7298A", "#A6761D","#D95F02" , "#E6AB02",  "#66A61E","#1B9E77", "#7570B3","#666666")
names(safe_colorblind_palette) <- levels(spe$precast_clusters)
for (pair in lr_pairs) {
  ligand <- pair[1]
  receptor <- pair[2]
  colname <- paste0("LR_", ligand, "_", receptor)
  print(colname)
  plot_list <- list()
  for (donor in sample_order) {
    spe_sub <- spe[ ,spe$sample_id == donor]
    spe_sub <- spe_sub[ ,!spe_sub$exclude_overlapping]
     spe_sub$coexpressing_spots <- factor(
    ifelse(spe_sub[[colname]] == "Coexpressing", "TRUE", "FALSE"),
    levels = c("TRUE", "FALSE")
    )
    plot_list[[donor]] <- make_escheR(spe_sub) |> 
                          add_fill(var = "coexpressing_spots") |> 
                          add_ground(var = "precast_clusters", stroke = 0.5) +
                          scale_color_manual(values = safe_colorblind_palette) +
                          scale_fill_manual(values = c("TRUE" = "#525252", "FALSE" = "white")) 
  }
  pdf(file.path(plot_dir, paste0(colname, "_spot_plots_escheR.pdf")), width = 12, height = 12)
  for (p in plot_list) print(p)
  dev.off()
}

for (pair in lr_pairs) {
  ligand <- pair[1]
  receptor <- pair[2]
  colname <- paste0("LR_", ligand, "_", receptor)
  spe$coexpressing_spots <- factor(
    ifelse(spe[[colname]] == "Coexpressing", "TRUE", "FALSE"),
    levels = c("TRUE", "FALSE")
    )
  plot_list <- list()
  for (donor in sample_order) {
    plot_list[[donor]] <- vis_clus(spe, sampleid = donor, clustervar = "coexpressing_spots", is_stitched = TRUE, 
    colors = c("red", "#CCCCCC40"))
  }
  spe$coexpressing_spots <- NULL
  pdf(file.path(plot_dir, paste0(colname, "_spot_plots.pdf")), width = 7, height = 7)
  for (p in plot_list) print(p)
  dev.off()
}

##########################################################################
# Add colocalization heatmap
##########################################################################
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
  "OPC" = "#0D9F72"
)
# Function to compute and plot co-localization matrix
# Function to compute and plot co-localization matrix
plot_LR_colocalization <- function(spe, ligand, receptor, cell_type_colors, min_expr = 1, outdir = plot_dir) {
  # 1. Check genes
  expr_mat <- assays(spe)$logcounts
  if (!ligand %in% rownames(expr_mat) || !receptor %in% rownames(expr_mat)) {
    warning(glue("Gene {ligand} or {receptor} not found in data."))
    return(NULL)
  }

  # 2. Expression-based spot categorization
  lig_expr <- expr_mat[ligand, ] > min_expr
  rec_expr <- expr_mat[receptor, ] > min_expr
  coexpr_spots <- lig_expr & rec_expr
  ctrl_spots <- !(coexpr_spots)

  if (sum(coexpr_spots) < 5) {
    warning("Fewer than 5 co-expressing spots.")
    return(NULL)
  }

  # 3. Cell type matrix
  weights_mat <- as.matrix(colData(spe)[, names(cell_type_colors)])
  weights_mat[is.na(weights_mat)] <- 0

  compute_adjacency <- function(mat, idx) {
    spot_mat <- mat[idx, , drop = FALSE]
    t(spot_mat) %*% spot_mat
  }

  normalize_adjacency <- function(adj) {
    adj2 <- adj
    diag(adj2) <- 0
    adj2[upper.tri(adj2)] <- 0
    adj / sum(adj2)
  }

  # 4. Compute adjacency matrices
  adj_LR <- compute_adjacency(weights_mat, coexpr_spots)
  adj_ctrl <- compute_adjacency(weights_mat, ctrl_spots)

  adj_LR_norm <- normalize_adjacency(adj_LR)
  adj_ctrl_norm <- normalize_adjacency(adj_ctrl)

  # 5. Ratio matrix
  ratio_mat <- adj_LR_norm / (adj_ctrl_norm + 1e-6)

  # 6. Heatmap annotations with correct mapping
  ha_row <- rowAnnotation(bar = rownames(ratio_mat), 
  col = list(bar = cell_type_colors), show_legend = FALSE, show_annotation_name = FALSE)

  ha_col <- HeatmapAnnotation(bar = colnames(ratio_mat), 
  col = list(bar = cell_type_colors), show_legend = FALSE, show_annotation_name = FALSE)

  # 7. Build heatmap
  ht <- Heatmap(
    ratio_mat,
    rect_gp = gpar(type = "none"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    row_names_side = "left",
    column_names_side = "bottom",
    col = inferno(100, direction = -1),
    name = "LR/Control",
    left_annotation = ha_row,
    bottom_annotation = ha_col,
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (i > j) {
        grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
      }else{
        grid.rect(x, y, w, h, gp = gpar(fill = "grey80", col = "grey80"))
      }
    }
  )

  # 8. Save to PDF
  pdf(file.path(outdir, paste0("colocalization_heatmap_", ligand, "_", receptor, ".pdf")), width = 5, height = 4.6)
  draw(ht, column_title = paste0("Colocalizing cell types in\n", ligand, " and ", receptor, " co-expressing spots"), column_title_gp = gpar(fontsize = 12))
  dev.off()
}

# Run for each LR pair
for (pair in lr_pairs) {
  print(pair)
  plot_LR_colocalization(spe, ligand = pair[1], receptor = pair[2], cell_type_colors)
}

plot_rctd_box_for_lr <- function(
  spe,
  lr_col      = "LR_PENK_OPRM1",
  celltypes   = c("DRD1_MSN_B","DRD1_MSN_D","DRD2_MSN_A"),
  cluster_col = "precast_clusters",
  restrict_domains = NULL,               # e.g., c("MSN 1","MSN 2","MSN 3","D1 islands")
  outline_color = "black",
  alpha = 0.9,
  notch = FALSE,
  varwidth = FALSE
) {
  stopifnot(lr_col %in% colnames(colData(spe)))
  stopifnot(all(celltypes %in% colnames(colData(spe))))

  # Parse ligand/receptor from lr_col for labels & colors
  parts <- strsplit(gsub("^LR_", "", lr_col), "_")[[1]]
  ligand   <- if (length(parts) >= 1) parts[1] else "Ligand"
  receptor <- if (length(parts) >= 2) parts[2] else "Receptor"

  # LR category order (y-axis, bottom->top after coord_flip)
  lr_levels <- c("Coexpressing", paste0(ligand, " only"), paste0(receptor, " only"), "None")
  y_labs <- setNames(
    c(paste0(ligand, " & ", receptor), paste0(ligand, " only"), paste0(receptor, " only"), "Neither"),
    lr_levels
  )

  # Colors for classes (coexpr / ligand-only / receptor-only); "Neither" uses na.value
  cat_labels <- c(paste0(ligand, " & ", receptor), paste0(ligand), paste0(receptor))
  cat_colors <- setNames(c("#d73027", "#91bfdb", "#4575b4"), cat_labels)

  library(dplyr); library(tidyr); library(ggplot2)

  # Data
  df <- as.data.frame(colData(spe)) |>
    dplyr::filter(is.null(restrict_domains) | .data[[cluster_col]] %in% restrict_domains) |>
    dplyr::filter(!is.na(.data[[lr_col]])) |>
    dplyr::mutate(
      Category = factor(.data[[lr_col]], levels = lr_levels),
      # mapping for fill colors:
      CatFill = dplyr::case_when(
        Category == "Coexpressing" ~ paste0(ligand, " & ", receptor),
        Category == paste0(ligand, " only") ~ paste0(ligand),
        Category == paste0(receptor, " only") ~ paste0(receptor),
        TRUE ~ NA_character_        # "None" -> NA -> uses na.value
      )
    ) |>
    dplyr::select(Category, CatFill, dplyr::all_of(celltypes)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(celltypes),
      names_to = "CellType",
      values_to = "RCTD_weight"
    ) |>
    dplyr::filter(is.finite(RCTD_weight))

  # Plot — horizontal boxplots; shared y, free x per facet
  ggplot(df, aes(x = Category, y = RCTD_weight, fill = CatFill)) +
    geom_boxplot(
      color = outline_color, alpha = alpha,
      notch = notch, varwidth = varwidth, outlier.size = 0.6, width = 0.7
    ) +
    facet_wrap(~ CellType, nrow = 1, scales = "free_x") +
    scale_x_discrete(labels = y_labs) +
    scale_fill_manual(values = cat_colors, breaks = cat_labels, na.value = "grey80", name = NULL) +
    labs(x = NULL, y = "RCTD weight (cell-type proportion)") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      axis.text.y = element_text(size = 10)
    ) +
    coord_flip()
}



p <- plot_rctd_box_for_lr(
  spe,
  lr_col = "LR_NLGN1_NRXN1", 
  celltypes   = c("DRD1_MSN_A", "DRD1_MSN_C", "OPC", "Astrocyte_A", "Excitatory", "DRD1_MSN_B", "Inh_B")
)

pdf(file.path(plot_dir, "RCTD_LR_NLGN1_NRXN1.pdf"), width = 12, height = 3)
print(p)
dev.off()