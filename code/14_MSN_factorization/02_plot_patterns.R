rm(list = ls())
library(MERINGUE)
library(sessioninfo)
library(here)
library(HDF5Array)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)
library(spatialNAcUtils)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)
library(tidyr)
library(purrr)
library(proxy)
library(vegan)
library(stats)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(scales)


# Set directories
plotDir <- here("plots", "14_MSN_factorization", "01_meringue")
resDir <- here("processed-data", "14_MSN_factorization", "01_meringue")
# Read in the processed SPE data
spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_dimRed_hdf5"
)
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

# Add clustering analysis from PRECAST
clusters_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters", "precast_clusters.csv")
spe[["spatial_domains"]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(clusters_path), by = 'key') |>
    pull(cluster) |>
    as.factor()
# Remove spots with no assigned clusters
spe <- spe[ ,!is.na(spe[["spatial_domains"]])]

# Format the spatialCoords to work with the gene expression data
rownames(spatialCoords(spe)) = colnames(spe)

# Subset to a particular sample if not transforming the row and column
#sample_id <- levels(spe$sample_id)[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]
sample_ids <- c("Br2743", "Br6432", "Br6423", "Br2720", "Br6471", "Br6522","Br8492", "Br8325", "Br8667", "Br3942")
spe_list <- lapply(sample_ids, function(sample_id){
  spe_subset <- spe[ ,spe$sample_id == sample_id]
  spe_subset <- spe_subset[ ,grepl("MSN", spe_subset[["spatial_domains"]])]
  cat(dim(spe_subset), "\n")
  spe_subset
})
names(spe_list) <- sample_ids

W_list <- lapply(spe_list, function(spe_subset){
  getSpatialNeighbors(spatialCoords(spe_subset), filterDist = 500)
})
I_list <- lapply(sample_ids, function(sample_id){
  readRDS(file.path(resDir, paste0("I_", sample_id, ".rds")))
})
names(I_list) <- sample_ids
results_filtered_list <- lapply(sample_ids, function(sample_id){
  readRDS(file.path(resDir, paste0("results_filtered_", sample_id, ".rds")))
})
names(results_filtered_list) <- sample_ids
scc_list <- lapply(sample_ids, function(sample_id){
  readRDS(file.path(resDir, paste0("scc_", sample_id, ".rds")))
})
names(scc_list) <- sample_ids
mat_list <- lapply(spe_list, function(spe_subset){
  mat <- logcounts(spe_subset)
  rownames(mat) <- make.names(rownames(mat), unique = TRUE)
  mat
})
safe_colorblind_palette <- c("#66A61E","#1B9E77", "#7570B3")
for(i in c(1:length(spe_list))){
  pdf(file.path(plotDir, paste0("Clusters_", sample_ids[i], ".pdf")))
  p <- spot_plot(
            spe_list[[i]],
            sample_id = sample_ids[i], var_name = "spatial_domains",
            is_discrete = TRUE, spatial = TRUE, colors = safe_colorblind_palette) + ggtitle(sample_ids[i])+
            #   Increase size of colored dots in legend
            guides(fill = guide_legend(override.aes = list(size = 5)))
  print(p)
  dev.off()
}

ggroup_list <- list()
for(i in c(1:length(spe_list))){
  spe_subset <- spe_list[[i]]
  mat <- mat_list[[i]]
  results.filter <- results_filtered_list[[i]]
  scc <- scc_list[[i]]
  pdf(file.path(plotDir, paste0("Patterns_", sample_ids[i], "_raw.pdf")))
  par(mfrow=c(2,2), mar=rep(2,4))
  ggroup <- groupSigSpatialPatterns(pos = spatialCoords(spe_subset), 
                                  mat = as.matrix(mat[results.filter,]), 
                                  scc = scc, 
                                  power = 1, 
                                  hclustMethod = 'ward.D', 
                                  deepSplit = 2,
                                  zlim=c(-1.5,1.5))
  dev.off()
  patterns <- do.call(cbind, ggroup[["prs"]])
  colnames(patterns) <- paste0("pattern_", c(1:dim(patterns)[2]))
  patterns <- patterns[match(rownames(colData(spe_subset)), rownames(patterns)), ]
  colData(spe_subset) <- cbind(colData(spe_subset), patterns)
  spe_list[[i]] <- spe_subset
  ggroup_list[[i]] <- ggroup
}

for(i in c(1:length(spe_list))){
  cat(sample_ids[i], "\n")
  nPatterns <- sum(grepl("pattern", colnames(colData(spe_list[[i]]))))
  plot_list <- list()
  for(j in c(1:nPatterns)){
    plot_list[[j]] <- vis_gene(spe_list[[i]], sampleid = sample_ids[i], geneid = paste0("pattern_", j), is_stitched = TRUE, cont_colors = viridisLite::magma(10, direction = -1)) + ggtitle(paste0("Pattern ", j))
  }
  pdf(file.path(plotDir, paste0("Patterns_", sample_ids[i], ".pdf")))
  for(j in c(1:nPatterns)){
   print(plot_list[[j]])
  }
  dev.off()
}

# --------------------------
# Step 1: Pattern Assigned Genes
# --------------------------
pattern_assigned_genes <- list()
for (i in seq_along(sample_ids)) {
  id <- sample_ids[i]
  groups <- ggroup_list[[i]]$groups
  for (pattern in levels(groups)) {
    key <- paste0(id, "_Pattern_", pattern)
    pattern_assigned_genes[[key]] <- names(groups[groups == pattern])
  }
}
saveRDS(pattern_assigned_genes, file.path(resDir, "all_pattern_assigned_genes.rds"))

# --------------------------
# Step 2: Domain Profiles
# --------------------------
domain_profiles <- list()
for (i in seq_along(sample_ids)) {
  id <- sample_ids[i]
  eigs <- ggroup_list[[i]]$prs
  spe_subset <- spe_list[[i]]
  spe_subset$spatial_domains <- factor(spe_subset$spatial_domains, levels = c("MSN 1", "MSN 2", "MSN 3"))
  for (j in seq_along(eigs)) {
    pattern_id <- paste0(id, "_Pattern_", names(eigs)[j])
    domain_profiles[[pattern_id]] <- tapply(eigs[[j]], spe_subset$spatial_domains, mean, na.rm = TRUE)
  }
}
saveRDS(domain_profiles, file.path(resDir, "domain_profiles.rds"))

# -------------------------------
# Step 1: Compute Similarity Matrices
# -------------------------------

compute_jaccard_matrix <- function(gene_list) {
  keys <- names(gene_list)
  n <- length(keys)
  mat <- matrix(NA, nrow = n, ncol = n)
  rownames(mat) <- colnames(mat) <- keys
  for (i in 1:n) {
    for (j in i:n) {
      inter <- length(intersect(gene_list[[i]], gene_list[[j]]))
      union <- length(union(gene_list[[i]], gene_list[[j]]))
      val <- ifelse(union == 0, 0, inter / union)
      mat[i, j] <- mat[j, i] <- val
    }
  }
  mat
}

get_corr_similarity_matrix <- function(list_input) {
  df <- bind_rows(lapply(list_input, function(x) tibble::enframe(x, name = "domain", value = "score")), .id = "pattern") |>
    pivot_wider(names_from = domain, values_from = score) |>
    column_to_rownames("pattern")
  df <- df[complete.cases(df), ]
  cor(t(df), use = "pairwise.complete.obs")
}

jaccard_pattern_genes <- compute_jaccard_matrix(pattern_assigned_genes)
domain_cor_matrix <- get_corr_similarity_matrix(domain_profiles)

# -------------------------------
# Step 2: Compute Consensus Similarity
# -------------------------------

pattern_keys <- intersect(rownames(jaccard_pattern_genes), rownames(domain_cor_matrix))
n <- length(pattern_keys)
consensus_sim <- matrix(NA, nrow = n, ncol = n)
rownames(consensus_sim) <- colnames(consensus_sim) <- pattern_keys

for (i in 1:n) {
  for (j in i:n) {
    x <- c(
      jaccard_pattern_genes[pattern_keys[i], pattern_keys[j]],
      domain_cor_matrix[pattern_keys[i], pattern_keys[j]]
    )
    val <- mean(x, na.rm = TRUE)
    consensus_sim[i, j] <- consensus_sim[j, i] <- val
  }
}

# -------------------------------
# Step 3: Cluster Patterns
# -------------------------------

hc <- hclust(as.dist(1 - consensus_sim), method = "ward.D2")
cluster_labels <- cutree(hc, k = 6)
pattern_clusters <- setNames(cluster_labels, pattern_keys)

# -------------------------------
# Step 4: Prepare Metadata + Colors
# -------------------------------

sample_ids <- sub("_Pattern_.*", "", pattern_keys)
annotation_df <- data.frame(
  Sample = sample_ids,
  Cluster = as.factor(pattern_clusters[pattern_keys])
)
rownames(annotation_df) <- pattern_keys

# Color palettes
sample_colors <- setNames(scales::hue_pal()(length(unique(annotation_df$Sample))), unique(annotation_df$Sample))
cluster_colors <- setNames(rainbow(length(unique(annotation_df$Cluster))), sort(unique(annotation_df$Cluster)))

ann_colors <- list(
  Sample = sample_colors,
  Cluster = cluster_colors
)

# -------------------------------
# Step 5: Plot Heatmaps
# -------------------------------

# Jaccard similarity heatmap
pdf(file.path(plotDir, "Jaccard_similarity_patterns.pdf"), width = 10, height = 10)
pheatmap(jaccard_pattern_genes[pattern_keys, pattern_keys],
         clustering_method = "ward.D2",
         color = colorRampPalette(c("white", "steelblue"))(100),
         annotation_row = data.frame(Sample = annotation_df$Sample),
         annotation_col = data.frame(Sample = annotation_df$Sample),
         annotation_colors = ann_colors,
         main = "Jaccard Index - Pattern Assigned Genes")
dev.off()

# Domain profile similarity heatmap
pdf(file.path(plotDir, "Domain_profile_correlation.pdf"), width = 10, height = 10)
pheatmap(domain_cor_matrix[pattern_keys, pattern_keys],
         clustering_method = "ward.D2",
         color = colorRampPalette(c("white", "darkgreen"))(100),
         annotation_row = data.frame(Sample = annotation_df$Sample),
         annotation_col = data.frame(Sample = annotation_df$Sample),
         annotation_colors = ann_colors,
         main = "Correlation - Domain Profiles")
dev.off()

# Consensus similarity heatmap
pdf(file.path(plotDir, "Consensus_similarity_patterns.pdf"), width = 10, height = 10)
pheatmap(consensus_sim,
         clustering_method = "ward.D2",
         color = colorRampPalette(c("white", "steelblue"))(100),
         annotation_row = data.frame(Cluster = annotation_df$Cluster),
         annotation_col = data.frame(Sample = annotation_df$Sample),
         annotation_colors = ann_colors,
         main = "Consensus Similarity - Patterns")
dev.off()
# --------------------------
# Grouped Boxplots: Pattern Values by MSN Domains
# --------------------------

# Output
pdf(file.path(plotDir, "PatternGrouped_Boxplot.pdf"), width = 5, height = 6)
for (sample_id in sample_ids) {
  spe_subset <- spe_list[[sample_id]]
  
  # Subset to relevant domains
  spe_sub <- spe_subset[, spe_subset$spatial_domains %in% c("MSN 1", "MSN 2", "MSN 3")]
  if (ncol(spe_sub) == 0) next
  
  # Extract pattern columns
  pattern_cols <- grep("^pattern_", colnames(colData(spe_sub)), value = TRUE)
  
  # Data frame for plotting
  df <- as.data.frame(colData(spe_sub)[, pattern_cols, drop = FALSE])
  df$Domain <- factor(as.character(spe_sub$spatial_domains), levels = c("MSN 1", "MSN 2", "MSN 3"))
  df$Spot <- colnames(spe_sub)
  
  # Pivot longer for ggplot
  df_long <- df %>%
    pivot_longer(cols = starts_with("pattern_"), names_to = "Pattern", values_to = "PatternValue")
  
  # Capitalize pattern labels for nicer x-axis display
  pattern_labels <- unique(df_long$Pattern)
  pretty_labels <- setNames(
    gsub("pattern_", "Pattern ", pattern_labels),
    pattern_labels
  )
  
  # Plot
  p <- ggplot(df_long, aes(x = Pattern, y = PatternValue, fill = Domain)) +
    geom_violin(position = position_dodge(0.75), width = 0.8, alpha = 0.6, trim = FALSE) +
    geom_boxplot(position = position_dodge(0.75), width = 0.2, outlier.shape = NA, color = "black") +
    scale_fill_manual(values = safe_colorblind_palette) +
    theme_classic(base_size = 12) +
    labs(title = paste("Sample:", sample_id),
         x = "Pattern",
         y = "Pattern Value") +
    scale_x_discrete(labels = pretty_labels) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")
  
  print(p)
}
dev.off()

# --------------------------
# Dotplot of Domain Activation per Pattern by Cluster
# --------------------------
dot_data <- domain_profiles %>%
  enframe("Pattern", "DomainMeans") %>% 
  unnest_wider(DomainMeans) %>% 
  rowwise() %>% 
  mutate(Sample = sub("_Pattern_.*", "", Pattern)) %>% 
  ungroup() %>%
  pivot_longer(cols = c("MSN 1", "MSN 2", "MSN 3"), 
               names_to = "Domain", values_to = "Mean") %>%
  mutate(
    Prop = map2_dbl(Pattern, Domain, function(p, d) {
      spe_subset <- spe_list[[ sub("_Pattern_.*", "", p) ]]
      v <- colData(spe_subset)[[ sub(".*_", "pattern_", p) ]]
      mean(v[spe_subset$spatial_domains == d] > 0, na.rm = TRUE)
    }),
    Cluster = factor(pattern_clusters[Pattern]),
    Mean = ifelse(is.na(Mean), 0, Mean)
  )

# --------------------------
# Step 2: Create mean and proportion matrices
# --------------------------

# Mean pattern values
mat_mean <- dot_data %>%
  select(Pattern, Domain, Mean) %>%
  pivot_wider(names_from = Domain, values_from = Mean) %>%
  column_to_rownames("Pattern") %>%
  as.matrix()

# Proportion of spots > 0
mat_prop <- dot_data %>%
  select(Pattern, Domain, Prop) %>%
  pivot_wider(names_from = Domain, values_from = Prop) %>%
  column_to_rownames("Pattern") %>%
  as.matrix()

# Order rows by cluster
pattern_order <- dot_data %>%
  distinct(Pattern, Cluster) %>%
  arrange(Cluster) %>%
  pull(Pattern)

mat_mean <- mat_mean[pattern_order, , drop = FALSE]
mat_prop <- mat_prop[pattern_order, , drop = FALSE]

# --------------------------
# Step 3: Color and size functions
# --------------------------

# Color scale for mean values
col_fun <- colorRamp2(
  c(min(mat_mean, na.rm = TRUE), 0, max(mat_mean, na.rm = TRUE)),
  c("purple", "white", "darkgreen")
)

# Dot size scaling function
scale_prop_to_size <- function(p) {
  rescale(p, to = c(1.5, 7), from = range(mat_prop, na.rm = TRUE))
}

# Dot size legend
legend_sizes <- scale_prop_to_size(c(0.25, 0.5, 0.75, 1.00)) * 2  # convert to diameter
dot_legend <- Legend(
  at = c(0.25, 0.5, 0.75, 1.00),
  type = "points",
  pch = 21,
  size = unit(legend_sizes, "pt"),  # diameter in pt
  legend_gp = gpar(fill = "grey60", col = "black"),
  title = "Proportion > 0"
)

# --------------------------
# Step 4: Row annotation for clusters
# --------------------------

cluster_factor <- factor(pattern_clusters[rownames(mat_mean)])
cluster_anno <- rowAnnotation(
  Cluster = cluster_factor,
  col = list(Cluster = setNames(rainbow(length(unique(cluster_factor))),
                                sort(unique(cluster_factor))))
)

# --------------------------
# Step 5: Heatmap plot
# --------------------------

# Create the Heatmap object
ht <- Heatmap(
  matrix = mat_mean,
  col = col_fun,
  name = "Mean Pattern\nValue",
  rect_gp = gpar(type = "none"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.circle(
      x = x,
      y = y,
      r = unit(scale_prop_to_size(mat_prop[i, j]), "pt"),
      gp = gpar(fill = col_fun(mat_mean[i, j]), col = "black")
    )
  },
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 12, fontface = "bold"),
  column_names_gp = gpar(fontsize = 14, fontface = "bold"),
  column_title = "MSN Domains",
  row_title = "Patterns (Grouped by Cluster)",
  left_annotation = cluster_anno,
  row_split = cluster_factor,
  row_gap = unit(2, "mm"),
  border = TRUE
)

# Output to PDF
pdf(file.path(plotDir, "Patterns_MSN_dotplot.pdf"), width = 6, height = 12)
draw(ht, annotation_legend_list = list(dot_legend))
dev.off()

# Filter to consensus clusters 1–4
consensus_subset <- pattern_clusters[pattern_clusters %in% 1:4]

# Step 1: Top genes per cluster
genes_by_cluster <- split(names(consensus_subset), consensus_subset)

gene_count_by_cluster <- lapply(genes_by_cluster, function(patterns) {
  gene_lists <- pattern_assigned_genes[patterns]
  data.frame(table(unlist(gene_lists)))
})

gene_count_by_cluster[[1]]$cluster <- "Cluster_1"
gene_count_by_cluster[[2]]$cluster <- "Cluster_2"
gene_count_by_cluster[[3]]$cluster <- "Cluster_3"
gene_count_by_cluster[[4]]$cluster <- "Cluster_4"

gene_count_by_cluster <- do.call(rbind, gene_count_by_cluster)
rownames(gene_count_by_cluster) <- c(1:dim(gene_count_by_cluster)[1])

colnames(gene_count_by_cluster) <- c("Gene", "Num_patterns", "MERINGUE_consensus")
write.csv(gene_count_by_cluster, file.path(resDir, "genes_per_meringue_consensus.csv"), row.names = FALSE, quote = FALSE)


# Sample order c("Br2743", "Br6432", "Br6423", "Br2720", "Br6471", "Br6522","Br8492", "Br8325", "Br8667", "Br3942")
pattern_df <- data.frame(
  Pattern = names(pattern_clusters),
  Cluster = as.integer(pattern_clusters)
)

# Extract sample ID
pattern_df$Sample <- sub("_Pattern_.*", "", pattern_df$Pattern)

# Focus only on clusters 1–4
pattern_df <- pattern_df[pattern_df$Cluster %in% 1:4, ]

# Count number of patterns per sample per cluster
pattern_counts <- pattern_df %>%
  group_by(Sample, Cluster) %>%
  tally()

# Keep only (Sample, Cluster) combinations with 1 pattern
unique_sample_cluster <- pattern_counts %>%
  filter(n == 1) %>%
  select(Sample, Cluster)

# Join back to get the pattern names
unique_patterns <- pattern_df %>%
  inner_join(unique_sample_cluster, by = c("Sample", "Cluster"))

unique_patterns <- rbind(unique_patterns, data.frame("Pattern" = "Br6432_Pattern_3", 
                                                     "Cluster" = 3, 
                                                     "Sample" = "Br6432"))

# Prepare a list to hold updated SPEs
spe_list_updated <- list()

# For each sample in spe_list
for (sample_id in names(spe_list)) {
  spe_sub <- spe_list[[sample_id]]
  col_data <- colData(spe_sub)
  
  # Get all patterns from this sample that are in unique_patterns
  sample_patterns <- unique_patterns %>% filter(Sample == sample_id)
  
  # Initialize empty columns for cluster_1 to cluster_4
  for (cl in 1:4) {
    col_data[[paste0("meringue_cluster_", cl)]] <- NA_real_
  }
  
  # For each selected pattern, assign its values to the corresponding cluster column
  for (i in seq_len(nrow(sample_patterns))) {
    pat_row <- sample_patterns[i, ]
    pattern_name <- pat_row$Pattern
    cluster_id <- pat_row$Cluster
    pattern_col <- sub(".*Pattern_", "pattern_", pattern_name)  # e.g., pattern_1
    
    # Only assign if the column exists (some samples may not have all patterns)
    if (pattern_col %in% colnames(col_data)) {
      cluster_col <- paste0("meringue_cluster_", cluster_id)
      col_data[[cluster_col]] <- col_data[[pattern_col]]
    }
  }
  
  # Update SPE with new colData
  colData(spe_sub) <- col_data
  spe_list_updated[[sample_id]] <- spe_sub
}

# Add MERINGUE consensus patterns to the original spe object
col_data <- lapply(spe_list_updated, function(spe_sub){
  colData(spe_sub)[ ,grepl("meringue", colnames(colData(spe_sub)))]
})
col_data <- do.call(rbind, col_data)
spe <- spe[ ,spe$spatial_domains %in% c("MSN 1", "MSN 2", "MSN 3")]
spe$spatial_domains <- factor(as.character(spe$spatial_domains), levels = c("MSN 1", "MSN 2", "MSN 3"))
col_data <- col_data[match(colnames(spe), rownames(col_data)), ]
colData(spe) <- cbind(colData(spe), col_data)


# ---- PARAMETERS ----
top_n <- 100  # Number of top genes per consensus pattern
pattern_cols <- grep("^meringue_cluster_", colnames(colData(spe)), value = TRUE)
expr_mat <- logcounts(spe)  # Can also use counts() or assays(spe)[["logcounts"]] if needed

for(pat in pattern_cols){
  spe_subset <- spe[ ,!is.na(spe[[pat]])]
  unique_samples <- as.character(unique(spe_subset$sample_id))
  plot_list <- list()
  for(j in unique_samples){
    plot_list[[j]] <- vis_gene(spe_subset, sampleid = j, geneid = pat, is_stitched = TRUE, cont_colors = viridisLite::magma(10, direction = -1)) + ggtitle(pat)
  }
  pdf(file.path(plotDir, paste0(pat, ".pdf")))
  for(j in unique_samples){
   print(plot_list[[j]])
  }
  dev.off()
}

# Ensure gene names are valid and unique
rownames(expr_mat) <- make.unique(rownames(expr_mat))

compute_cor_matrix <- function(expr, patterns) {
  top_genes <- list()
  cor_mat <- matrix(NA, nrow = nrow(expr), ncol = length(patterns))
  rownames(cor_mat) <- rownames(expr)
  colnames(cor_mat) <- patterns

  for (pat in patterns) {
    print(pat)
    pat_values <- colData(spe)[[pat]]
    non_na_idx <- which(!is.na(pat_values))
    expr_sub <- expr[, non_na_idx, drop = FALSE]
    pat_sub <- pat_values[non_na_idx]
    expr_sub <- as.matrix(expr_sub)
    cors <- apply(expr_sub, 1, function(x) cor(x, pat_sub, method = "pearson"))
    cor_mat[, pat] <- cors}

  cor_mat
}

cor_mat <- compute_cor_matrix(expr_mat, pattern_cols)

ranked_genes <- list()
for(pat in pattern_cols){
  cors <- cor_mat[ ,pat]
  ranked_genes[[pat]] <- names(sort(cors, decreasing = TRUE))
}

top_genes <- list()
for(pat in pattern_cols){
  top_genes[[pat]] <- ranked_genes[[pat]][1:top_n]
}

select_genes <- c("DRD1", "DRD2", "GAD1", "GAD2",
                  "PENK", "ARPP21", "PTPN5", "PDYN", "GNAL", "SLC32A1", "PRKCG", 
                  "MBP", "MOBP", "CLDND1", "OPALIN", "OLIG1", "GFAP", 
                  "PEG10", "VSNL1", "CALY", "NNAT", "NDN", "CARTPT",
                  "ANO3", "RASD2", "PDE1B", "ADORA2A", "PDE10A", "CALB1", "TAC1", "PPP1R1A", "ADCY5", 
                  "CPNE5", "CRYM","TSPAN7","SST", "NPY")

genes_available <- intersect(select_genes, rownames(cor_mat))
heat_mat <- cor_mat[genes_available, ,drop = FALSE]

# --- Heatmap color function ---
col_fun <- colorRamp2(c(-0.6, 0, 0.7), c("darkblue", "white", "darkred"))

# Define original rownames
original_order <- colnames(heat_mat)

# Create new labels and desired order
col_labels <- c(
  "meringue_cluster_4" = "MERINGUE Consensus Pattern 4",
  "meringue_cluster_1" = "MERINGUE Consensus Pattern 1",
  "meringue_cluster_3" = "MERINGUE Consensus Pattern 3",
  "meringue_cluster_2" = "MERINGUE Consensus Pattern 2"
)

# Reorder heat matrix rows
desired_order <- names(col_labels)
heat_mat <- heat_mat[ ,desired_order, drop = FALSE]

# Apply new names
colnames(heat_mat) <- col_labels[colnames(heat_mat)]

# --- Plot heatmap ---
pdf(file.path(plotDir, "meringue_consensus_selected_genes_heatmap.pdf"), width = 15, height = 4)
Heatmap(
  t(heat_mat),
  name = "Pearson\nCorrelation",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title = "Selected Genes",
  row_title = "MERINGUE Consensus Patterns",
  rect_gp = gpar(col = "black", lwd = 0.5),
  heatmap_legend_param = list(title = "Correlation", at = c(-0.6, 0, 0.7))
)
dev.off()