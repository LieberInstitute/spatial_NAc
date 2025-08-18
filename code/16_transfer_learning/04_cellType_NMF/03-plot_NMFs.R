library(tidyverse)
library(RcppML)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(here)
library(sessioninfo)
library(scater)
library(scran)
library(BiocParallel)
library(BiocSingular)
library(spatialNAcUtils)
library(jaffelab)
library(projectR)
library(scater)
library(scran)
library(dittoSeq)
library(escheR)
library(getopt)
library(Seurat)

opt <- list()
opt$data <- "rat_case_control_cocaine_acute"
opt$nFactors <- 30
opt$select_HVGs <- FALSE

res_dir <- here::here("processed-data", "16_transfer_learning", "04_cellType_NMF")
plot_dir <- here::here("plots", "16_transfer_learning", "04_cellType_NMF")
res_dir <- paste0(res_dir, "/", opt$data)

if(opt$select_HVGs){
    spe <- readRDS(file.path(res_dir, paste0("spe_NMF_", opt$nFactors ,"_HVGs.rds")))
}else{
    spe <- readRDS(file.path(res_dir, paste0("spe_NMF_", opt$nFactors ,".rds")))
}

plot_dir <- paste0(plot_dir, "/", opt$data, "/nmf_", opt$nFactors, "/")


###################################################################################################
# Plotting from below here
theme_Publication <- function(base_size=14, base_family="sans") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line.x = element_line(colour="black"),
               axis.line.y = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.box = "vetical",
               legend.key.size= unit(0.5, "cm"),
               #legend.margin = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
       ))
      
}

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#f87f01","#7fc97f","#ef3b2c","#feca01","#a6cee3","#fb9a99","#984ea3","#8C591D")), ...)
      
}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#f87f01","#7fc97f","#ef3b2c","#feca01","#a6cee3","#fb9a99","#984ea3","#8C591D")), ...)}
     
nmf_factors <- paste0("nmf", c(1:opt$nFactors))
nSpot_nonzero_coeff <- c()
for(i in c(1:length(nmf_factors))){
    nSpot_nonzero_coeff[i] <- sum(spe[[nmf_factors[i]]] > 0)
}
nSpot_nonzero_coeff_by_sample <- lapply(c(1:opt$nFactors), function(i){
    table(spe$sample_id[spe[[nmf_factors[i]]] > 0])
})
nSpot_nonzero_coeff_by_sample <- do.call(rbind, nSpot_nonzero_coeff_by_sample)
nSpot_nonzero_coeff_by_sample <- data.frame(nSpot_nonzero_coeff_by_sample)
nSpot_nonzero_coeff_by_sample$NMF <- c(1:opt$nFactors)
nSpot_nonzero_coeff_by_sample <- reshape2::melt(nSpot_nonzero_coeff_by_sample, id.vars = "NMF")

res_df <- data.frame("NMF" = nmf_factors, "num_spots" = nSpot_nonzero_coeff)
res_df$NMF <- factor(res_df$NMF, levels = nmf_factors)

res_df_final <- res_df[res_df$num_spots > 200, ]

spe_Br2743 <- mirrorObject(spe, sample_id = "Br2743", image_id = "lowres", axis = "v")
spe_Br8492 <- mirrorObject(spe, sample_id = "Br8492", image_id = "lowres", axis = "v")
spe_Br8325 <- mirrorObject(spe, sample_id = "Br8325", image_id = "lowres", axis = "v")
spe_Br3942 <- mirrorObject(spe, sample_id = "Br3942", image_id = "lowres", axis = "v")

spe <- spe[ ,!spe$sample_id %in% c("Br2743", "Br8492", "Br8325", "Br3942")]
spe <- cbind(spe, spe_Br2743)
spe <- cbind(spe, spe_Br8492)
spe <- cbind(spe, spe_Br8325)
spe <- cbind(spe, spe_Br3942)
sample_order <- c("Br2743", "Br6432", "Br6423", "Br2720", "Br6471", "Br6522","Br8492", "Br8325", "Br8667", "Br3942")

# Add clustering data
safe_colorblind_palette <- c("#E7298A", "#A6761D","#D95F02" , "#E6AB02",  "#66A61E","#1B9E77", "#7570B3","#666666")
clustDir <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters")
clusters_df <- read.csv(file.path(clustDir, "precast_clusters.csv"))
spe <- spe[ ,colnames(spe) %in% clusters_df$key]
clusters_df <- clusters_df[match(colnames(spe), clusters_df$key), ]
spe$PRECAST_clusters <- clusters_df$cluster
spe$PRECAST_clusters <- factor(spe$PRECAST_clusters, 
levels = c("D1 islands", "Endothelial/Ependymal", "Excitatory", "Inhibitory", "MSN 1", "MSN 2", "MSN 3", "WM"))

if(opt$data == "rat_case_control_cocaine_acute"){
    if(opt$nFactors == 30){
        select_nmfs <- c("nmf15", "nmf28")
    }
    if(opt$nFactors == 50){
        select_nmfs <- c("nmf17", "nmf42")
    }
}
if(opt$data == "rat_case_control_morphine_acute"){
    if(opt$nFactors == 30){
        select_nmfs <- c("nmf12", "nmf3")
    }
    if(opt$nFactors == 50){
        select_nmfs <- c("nmf30", "nmf7")
    }
}
if(opt$data == "rat_case_control_morphine_repeated"){
    if(opt$nFactors == 30){
        select_nmfs <- c("nmf12", "nmf3", "nmf23")
    }
    if(opt$nFactors == 50){
        select_nmfs <- c("nmf12", "nmf5", "nmf25", "nmf16")
    }
}

for(n in select_nmfs){
    cat(n, "\n")
     plot_donor_list <- lapply(sample_order, function(isample){
        cat(isample, "\n")
         vis_gene(spe, sampleid = isample, geneid = n, is_stitched = TRUE, minCount =0) + ggtitle(isample)
    })
    pdf(file.path(plot_dir, paste0(n, ".pdf")))
    print(plot_donor_list)
    dev.off()
}

plot_nmf_by_cluster <- function(df, nmf_column, group_column) {
  # Ensure the column names are unquoted symbols
  nmf_sym <- rlang::sym(nmf_column)
  group_sym <- rlang::sym(group_column)
  
  # Prepare data
  df_subset <- df %>%
    select(!!group_sym, !!nmf_sym) %>%
    filter(!is.na(!!group_sym))

  # Convert group column to factor
  df_subset[[group_column]] <- as.factor(df_subset[[group_column]])

  # Generate plot
  p <- ggplot(df_subset, aes_string(x = group_column, y = nmf_column, fill = group_column)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.7) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste(nmf_column, "by", group_column),
      x = group_column,
      y = nmf_column
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) + scale_fill_manual(values = safe_colorblind_palette)

    return(p)
}

df <- data.frame(colData(spe))
for(nmf in select_nmfs){
    p <- plot_nmf_by_cluster(df, nmf, "PRECAST_clusters")
    pdf(file.path(plot_dir, paste0(nmf,"_violin.pdf")), width = 8, height = 3)
    print(p)
    dev.off()
}

### Make dotplots
# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Build data frame for all combinations of selected NMFs and spatial domains
dotplot_df <- expand.grid(
  nmf = select_nmfs,
  cluster = levels(spe$PRECAST_clusters)
)

# Calculate proportion of non-zero values and mean projection for each NMF × cluster
dotplot_summary <- dotplot_df %>%
  rowwise() %>%
  mutate(
    values = list(spe[[nmf]][spe$PRECAST_clusters == cluster]),
    nonzero = sum(values[[1]] > 0, na.rm = TRUE),
    total = length(values[[1]]),
    proportion = ifelse(total > 0, nonzero / total, NA_real_),
    mean_value = ifelse(total > 0, mean(values[[1]], na.rm = TRUE), NA_real_)
  ) %>%
  ungroup() %>%
  select(nmf, cluster, proportion, mean_value)

  # Drop NAs
dotplot_summary <- na.omit(dotplot_summary)

# Make factors ordered for clean plotting
dotplot_summary$nmf <- factor(dotplot_summary$nmf, levels = select_nmfs)
dotplot_summary$cluster <- factor(dotplot_summary$cluster, levels = levels(spe$PRECAST_clusters))

# Plot dotplot

p_dot <- ggplot(dotplot_summary, aes(x = nmf, y = cluster)) +
  geom_point(aes(size = proportion, color = mean_value)) +
  scale_color_viridis_c(option = "inferno", direction = -1) +  # reversed inferno
  scale_size(range = c(1, 6)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "NMF Factor",
    y = "Spatial Domain",
    color = "Mean\nProjection",
    size = "Proportion\n> 0"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # box around plot
    panel.background = element_blank()
  )


# Save plot
pdf(file.path(plot_dir, "selected_nmf_by_cluster_dotplot.pdf"), width = 5, height = 5)
print(p_dot)
dev.off()

# Make correlation heatmap
compute_visium_gene_correlations <- function(spe, nmf_factors, assay_name = "logcounts") {
  expr <- assay(spe, assay_name)
  coldata <- colData(spe)

  cor_mat <- sapply(nmf_factors, function(fac) {
    vec <- coldata[[fac]]
    apply(expr, 1, function(g) cor(g, vec, method = "pearson"))
  })
   rownames(cor_mat) <- rownames(expr)
  return(cor_mat)
}

cor_mat <- compute_visium_gene_correlations(spe, select_nmfs)
rownames(cor_mat) <- make.unique(rowData(spe)$gene_name)


plot_visium_gene_correlation_heatmap <- function(cor_mat, select_nmfs, top_n = 20, output_dir, file_name) {
  library(dplyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)

  # 1. Select top_n genes per factor
  top_genes_list <- lapply(select_nmfs, function(fac) {
    if (!fac %in% colnames(cor_mat)) {
      stop(paste("Factor", fac, "not found in cor_mat"))
    }
    sort_genes <- sort(cor_mat[, fac], decreasing = TRUE)
    names(head(sort_genes, top_n))
  })
  names(top_genes_list) <- select_nmfs

  all_top_genes <- unique(unlist(top_genes_list))
  cor_mat_sub <- cor_mat[all_top_genes, select_nmfs, drop = FALSE]
  #cor_mat_sub <- t(cor_mat_sub)  # Rows = factors, Columns = genes

  # 2. Build gene-to-factor mapping
  gene_to_factor <- do.call(rbind, lapply(names(top_genes_list), function(fac) {
    genes <- top_genes_list[[fac]]
    if (length(genes) == 0) return(NULL)
    data.frame(factor = fac, gene = genes)
  }))

  if (is.null(gene_to_factor) || nrow(gene_to_factor) == 0) {
    stop("No top genes found for any NMF factor.")
  }

  gene_factor_labels <- gene_to_factor %>%
    group_by(gene) %>%
    summarise(Factor = paste(sort(unique(factor)), collapse = ", ")) %>%
    column_to_rownames("gene")

  # Keep only observed combinations (fix: no unused combinations)
  observed_levels <- c("nmf12, nmf23", "nmf12, nmf3", "nmf23", "nmf12", "nmf3")
  gene_factor_labels$Factor <- factor(gene_factor_labels$Factor, levels = observed_levels)
  gene_factor_labels <- gene_factor_labels[!is.na(gene_factor_labels$Factor), , drop = FALSE]
  cor_mat_sub <- cor_mat_sub[rownames(gene_factor_labels), , drop = FALSE]

  # 4. Define annotation colors
  n_groups <- length(levels(gene_factor_labels$Factor))
  palette_colors <- colorRampPalette(brewer.pal(8, "Set3"))(n_groups)
  factor_colors <- setNames(palette_colors, levels(gene_factor_labels$Factor))

  row_ha <- rowAnnotation(
    Factor = gene_factor_labels$Factor,
    col = list(Factor = factor_colors),
    show_annotation_name = TRUE
  )

  # 5. Color mapping for correlation values
  cor_range <- range(cor_mat_sub, na.rm = TRUE)
  breaks <- pretty(cor_range, n = 7)
  cor_colors <- circlize::colorRamp2(breaks, colorRampPalette(c("blue", "white", "red"))(length(breaks)))

  # 6. Adjust font size dynamically
  n_genes <- nrow(cor_mat_sub)
  label_fontsize <- 10

  # 7. Draw heatmap
  ht <- Heatmap(
    cor_mat_sub,
    name = "Correlation",
    col = cor_colors,
    left_annotation = row_ha,
    cluster_columns = TRUE,
    cluster_rows = FALSE,
    row_split = gene_factor_labels$Factor,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = label_fontsize),
    row_title = "Genes",
    column_title = "Factors",
    heatmap_legend_param = list(title = "Correlation")
  )

  # 8. Export to PDF
  pdf(file.path(output_dir, file_name), width = 5, height = 7)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       column_title = "Gene Correlation with NMF Projections",
       column_title_gp = gpar(fontsize = 14, fontface = "bold"))
  dev.off()
}

plot_visium_gene_correlation_heatmap(
  cor_mat = cor_mat,
  select_nmfs = select_nmfs,
  top_n = 20,
  output_dir = plot_dir,
  file_name = "visium_gene_projection_heatmap_MSN.pdf"
)

# Also look at genes by top loadings
res_dir <- here::here("processed-data", "16_transfer_learning", "04_cellType_NMF")
res_dir <- paste0(res_dir, "/", opt$data)
if(opt$select_HVGs){
    x <- readRDS(file.path(res_dir, paste0("nmf_results_HVGs_", opt$nFactors, "factors.rds")))
}else{
    x <- readRDS(file.path(res_dir, paste0("nmf_results_", opt$nFactors, "factors.rds")))
}
patterns <- t(x@h)
colnames(patterns) <- paste("NMF", 1:dim(patterns)[2], sep = "_")
# extract loadings
loadings <- x@w

loadings_sub <- loadings[ ,select_nmfs]

loadings_long <- loadings_sub %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = starts_with("nmf"), names_to = "factor", values_to = "loading")

# Filter top genes per factor
top_loadings <- loadings_long %>%
  group_by(factor) %>%
  top_n(20, wt = loading) %>%
  ungroup()

# Optional: reorder genes within each facet
top_loadings <- top_loadings %>%
  group_by(factor) %>%
  mutate(gene = fct_reorder(gene, loading)) %>%
  ungroup()

# Plot
pdf(file.path(plot_dir, "top_loadings_barplot.pdf"), width = 12, height = 6)
ggplot(top_loadings, aes(x = gene, y = loading, fill = factor)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(~factor, scales = "free") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
    panel.spacing = unit(1, "lines")
  ) + 
  labs(
    title = paste("Top", top_n, "Gene Loadings per NMF Factor"),
    x = "Gene",
    y = "Loading"
  )
dev.off()