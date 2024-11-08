library(here)
library(PRECAST)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(spatialNAcUtils)
library(purrr)
library(ggpubr)
library(ggsci)
library(dittoSeq)
library(getopt)
library(pheatmap)

spec <- matrix(
    c("nnSVG_type", "n", 1, "logical", "Use nnSVGs identified by controlling for PRECAST k = 2 clusters?", 
      "use_random_start", "r", 1, "logical", "Were random starts run for this version?"),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

if(opt$use_random_start){
    random_start <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}

if(opt$nnSVG_type){
    if(opt$use_random_start){
        out_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", paste0("random_start_", random_start), "PRECAST_k%s.csv")
        plot_dir <- here("plots", "07_spatial_domains", "01_precast", "nnSVG_precast", paste0("random_start_", random_start))
    }else{
        out_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "PRECAST_k%s.csv")
        plot_dir <- here("plots", "07_spatial_domains", "01_precast", "nnSVG_precast")
    }
  
}else{
     if(opt$use_random_start){
        out_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_default", paste0("random_start_", random_start), "PRECAST_k%s.csv")
        plot_dir <- here("plots", "07_spatial_domains", "01_precast", "nnSVG_default", paste0("random_start_", random_start))
    }else{
        out_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_default", "PRECAST_k%s.csv")
        plot_dir <- here("plots", "07_spatial_domains", "01_precast", "nnSVG_default")
    }
}

spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_dimRed_hdf5"
)


#   List of genes provided by Svitlana, named by the subregion they're markers
#   for
subregion_genes <- list(
    "shell" = c(
        "TAC3", "TAC1", "CALB2", "GRIA4", "CPNE4", "GREB1L", "ARHGAP6", "OPRK1"
    ),
    "core" = toupper(
        c(
            "CALB1", "HPCAL4", "ZDBF2", "LAMP5", "Scn4b", "Rgs9", "Ppp3ca",
            "Oprd1", "Adora2a", "Pde1a", "Peg10", "Dlk1", "Htr2c"
        )
    ),
    "white_matter" = c("MBP", "GFAP", "PLP1", "AQP4"),
    "D1_islands" = c("OPRM1", "DRD1", "RXFP1", "CPNE4", "FOXP2")
)

dir.create(plot_dir, showWarnings = FALSE)

spe <- loadHDF5SummarizedExperiment(spe_dir)
stopifnot(all(unlist(subregion_genes) %in% rowData(spe)$gene_name))

################################################################################
#   Import PRECAST results and append to colData
################################################################################

result_list <- list()
for (K in 3:12) {
    result_list[[K]] <- sprintf(out_path, K) |>
        read.csv() |>
        as_tibble() |>
        select(c(key, cluster)) |>
        mutate(k = K)
}
precast_results <- do.call(rbind, result_list) |>
    pivot_wider(
        values_from = cluster, names_from = k, names_prefix = "precast_k"
    )

precast_results <- data.frame(precast_results)
for(i in c(2:dim(precast_results)[2])){
    precast_results[ ,i] <- factor(as.character(precast_results[ ,i]), levels = c(1:(i+1)))
}

temp <- colnames(spe)
colData(spe) <- colData(spe) |>
    as_tibble() |>
    left_join(precast_results, by = "key") |>
    DataFrame()
colnames(spe) <- temp

################################################################################
#   Simply plot cluster assignments for k=2 through k=28
################################################################################

for (k in 3:12) {
    plot_list <- list()
    for (donor in unique(spe$sample_id)) {
        plot_list[[donor]] <- spot_plot(
                spe,
                sample_id = donor, var_name = paste0("precast_k", k),
                is_discrete = TRUE, spatial = TRUE
            ) +
            #   Increase size of colored dots in legend
            guides(fill = guide_legend(override.aes = list(size = 5)))
    }

    pdf(file.path(plot_dir, sprintf("k%s_all_samples.pdf", k)))
    print(plot_list)
    dev.off()
}

# For each value of k, determine the number of spots that belong to each cluster from each donor
p_list <- list()
index <- 1

for(k in 3:12){
    cluster_assignments <- data.frame(spe@colData)
    group_ids <- c("sample_id", paste0("precast_k", k))
    k_summary <- cluster_assignments %>%  group_by(across(all_of(group_ids))) %>% summarize(n = n())
    colnames(k_summary) <- c("sample_id", "cluster", "n")
    k_summary$cluster <- as.factor(k_summary$cluster)
    p_list[[index]] <- ggplot(k_summary, aes(x = cluster, y = n, fill = sample_id)) + geom_bar(position="stack", stat="identity") + 
    theme_classic(base_size = 22) + scale_fill_npg() + xlab("Cluster") + ylab("Number of spots") + guides(fill=guide_legend(title="Individual"))
    index <- index + 1
}

pdf(file.path(plot_dir, "cluster_spots_by_sample.pdf"), width = 15, height = 5)
for(k in 1:length(p_list)){
    print(p_list[[k]])
}
dev.off()

# Add in clustering for k = 2
WM_cluster_path <- here("processed-data", "07_spatial_domains","01_precast", "00_pre_clustering", "PRECAST_k2.csv")
WM_results <- read.csv(WM_cluster_path) |>
        as_tibble()
WM_cluster <- colData(spe) |>
        as_tibble() |>
        left_join(WM_results, by = "key") |>
        pull(cluster)
spe$precast_k2 <- WM_cluster
spe$WM_cluster <- WM_cluster
spe$WM_cluster[WM_cluster == 1] <- "Gray matter"
spe$WM_cluster[WM_cluster == 2] <- "White matter"
spe$WM_cluster <- as.factor(spe$WM_cluster)

spe <- spe[, !is.na(spe$WM_cluster)]

# Visualize proportion of white matter spots in each cluster
p_list <- list()
index <- 1
for(k in 3:12){
    prop_WM <- list()
    for(j in 1:k){
        ind <- spe@colData[ ,paste0("precast_k", k)] == j
        prop_WM[[j]] <- table(spe$WM_cluster[ind])
    }
    prop_WM <- data.frame(do.call(rbind, prop_WM))
    prop_WM$Cluster <- c(1:k)
    prop_WM <- reshape2::melt(prop_WM, id.vars = "Cluster")
    colnames(prop_WM) <- c("Cluster", "Region", "nSpots")
    prop_WM$Cluster <- factor(prop_WM$Cluster, levels = c(1:k))
    levels(prop_WM$Region) <- gsub("\\.", " ", levels(prop_WM$Region))
    p_list[[index]] <- ggplot(prop_WM, aes(x = Cluster, y = nSpots, fill = Region)) + geom_bar(stat = "identity", position = "stack") + 
    theme_classic(base_size = 22) + ylab("Number of spots") 
    index <- index + 1
}
pdf(file.path(plot_dir, "prop_WM_in_clusters.pdf"), width = 15, height = 5)
for(k in 1:length(p_list)){
    print(p_list[[k]])
}
dev.off()

# Plot the distribution of UMI counts in each cluster
p_list <- list()
index <- 1
for(k in c(3:12)){
  p_list[[index]] <- ggplot(data = as.data.frame(colData(spe)),
       aes_string(x = paste0("precast_k", k), y = "sum_umi", fill = paste0("precast_k", k))) + 
  geom_violin() + xlab("Cluster") + ylab("Library size") +
  theme_classic(base_size = 22) + scale_fill_igv() + ylim(c(0, 18000)) + ggtitle(paste0("K = ", k)) +
  theme(legend.position = "none")
  index <- index + 1
}
pdf(file.path(plot_dir, "Library_size_distributions_by_cluster.pdf"), width = 20, height = 5)
for(k in 1:length(p_list)){
    print(p_list[[k]])
}
dev.off()

# Distributions of number of genes
p_list <- list()
index <- 1
for(k in c(3:12)){
  p_list[[index]] <- ggplot(data = as.data.frame(colData(spe)),
       aes_string(x = paste0("precast_k", k), y = "sum_gene", fill = paste0("precast_k", k))) + 
  geom_violin() + xlab("Cluster") + ylab("Number of genes") +
  theme_classic(base_size = 22) + scale_fill_igv() + ylim(c(0, 10000)) + ggtitle(paste0("K = ", k)) +
  theme(legend.position = "none")
  index <- index + 1
}
pdf(file.path(plot_dir, "Detected_gene_distributions_by_cluster.pdf"), width = 20, height = 5)
for(k in 1:length(p_list)){
    print(p_list[[k]])
}
dev.off()

# Distributions of edge distance
p_list <- list()
index <- 1
for(k in c(3:12)){
  p_list[[index]] <- ggplot(data = as.data.frame(colData(spe)),
       aes_string(x = paste0("precast_k", k), y = "edge_distance", fill = paste0("precast_k", k))) + 
  geom_violin() + xlab("Cluster") + ylab("Edge distance") +
  theme_classic(base_size = 22) + scale_fill_igv() + ylim(c(0, 50)) + ggtitle(paste0("K = ", k)) +
  theme(legend.position = "none")
  index <- index + 1
}
pdf(file.path(plot_dir, "Edge_distance_distributions_by_cluster.pdf"), width = 20, height = 7)
for(k in 1:length(p_list)){
    print(p_list[[k]])
}
dev.off()

# Distribution of percent mito
p_list <- list()
index <- 1
for(k in c(3:12)){
  p_list[[index]] <- ggplot(data = as.data.frame(colData(spe)),
       aes_string(x = paste0("precast_k", k), y = "expr_chrM_ratio", fill = paste0("precast_k", k))) + 
  geom_violin() + xlab("Cluster") + ylab("Ratio of mitochondrial expression") +
  theme_classic(base_size = 22) + scale_fill_igv() + ylim(c(0, 0.6)) + ggtitle(paste0("K = ", k)) +
  theme(legend.position = "none")
  index <- index + 1
}
pdf(file.path(plot_dir, "percent_mito_distributions_by_cluster.pdf"), width = 20, height = 7)
for(k in 1:length(p_list)){
    print(p_list[[k]])
}
dev.off()

# Distribution of number of cells by cluster
p_list <- list()
index <- 1
for(k in c(3:12)){
  p_list[[index]] <- ggplot(data = as.data.frame(colData(spe)),
       aes_string(x = paste0("precast_k", k), y = "Nmask_dark_blue", fill = paste0("precast_k", k))) + 
  geom_violin() + xlab("Cluster") + ylab("Number of cells ") +
  theme_classic(base_size = 22) + scale_fill_igv() + ylim(c(0, 10)) + ggtitle(paste0("K = ", k)) +
  theme(legend.position = "none")
  index <- index + 1
}
pdf(file.path(plot_dir, "nCells_distributions_by_cluster.pdf"), width = 20, height = 7)
for(k in 1:length(p_list)){
    print(p_list[[k]])
}
dev.off()

# Check the distribution of shell genes in clusters
rownames(spe) <- rowData(spe)$gene_name

p_list <- list()
index <- 1
for(k in c(3:12)){
    cat(k, "\n")
    p_list[[index]] <- dittoPlot(spe, subregion_genes$shell, group.by = paste0("precast_k", k), plots = "vlnplot", assay = "logcounts", vlnplot.lineweight = 0.3)
    index <- index + 1
}
pdf(file.path(plot_dir, "shell_marker_expression.pdf"), width = 20, height = 10)
for(k in 1:length(p_list)){
    print(p_list[[k]])
}
dev.off()

p_list <- list()
index <- 1
for(k in c(3:12)){
    cat(k, "\n")
    p_list[[index]] <- dittoPlot(spe, subregion_genes$core, group.by = paste0("precast_k", k), plots = "vlnplot", assay = "logcounts", vlnplot.lineweight = 0.3)
    index <- index + 1
}
pdf(file.path(plot_dir, "core_marker_expression.pdf"), width = 20, height = 10)
for(k in 1:length(p_list)){
    print(p_list[[k]])
}
dev.off()

p_list <- list()
index <- 1
for(k in c(3:12)){
    cat(k, "\n")
    p_list[[index]] <- dittoPlot(spe, subregion_genes$white_matter, group.by = paste0("precast_k", k), plots = "vlnplot", assay = "logcounts", vlnplot.lineweight = 0.3)
    index <- index + 1
}
pdf(file.path(plot_dir, "WM_marker_expression.pdf"), width = 20, height = 10)
for(k in 1:length(p_list)){
    print(p_list[[k]])
}
dev.off()

p_list <- list()
index <- 1
for(k in c(3:12)){
    cat(k, "\n")
    p_list[[index]] <- dittoPlot(spe, subregion_genes$D1_islands, group.by = paste0("precast_k", k), plots = "vlnplot", assay = "logcounts", vlnplot.lineweight = 0.3)
    index <- index + 1
}
pdf(file.path(plot_dir, "D1_islands_marker_expression.pdf"), width = 20, height = 10)
for(k in 1:length(p_list)){
    print(p_list[[k]])
}
dev.off()


## Heatmap to visualize expression in clusters
p_list <- list()
index <- 1
for(k in c(3:12)){
    cat(k, "\n")
    p_list[[index]] <- dittoHeatmap(spe, subregion_genes$shell, group.by = paste0("precast_k", k), 
     order.by = paste0("precast_k", k), scale = "none", heatmap.colors = colorRampPalette(c("white", "red"))(50))
    index <- index + 1
}
pdf(file.path(plot_dir, "shell_marker_expression_heatmap.pdf"), width = 15, height = 10)
for(k in 1:length(p_list)){
    print(p_list[[k]])
}
dev.off()

p_list <- list()
index <- 1
for(k in c(3:12)){
    cat(k, "\n")
    p_list[[index]] <- dittoHeatmap(spe, subregion_genes$core, group.by = paste0("precast_k", k), 
     order.by = paste0("precast_k", k), scale = "none", heatmap.colors = colorRampPalette(c("white", "red"))(50))
    index <- index + 1
}
pdf(file.path(plot_dir, "core_marker_expression_heatmap.pdf"), width = 20, height = 7)
for(k in 1:length(p_list)){
    print(p_list[[k]])
}
dev.off()

p_list <- list()
index <- 1
for(k in c(3:12)){
    cat(k, "\n")
    p_list[[index]] <- dittoHeatmap(spe, subregion_genes$white_matter, group.by = paste0("precast_k", k), 
     order.by = paste0("precast_k", k), scale = "none", heatmap.colors = colorRampPalette(c("white", "red"))(50))
    index <- index + 1
}
pdf(file.path(plot_dir, "WM_marker_expression_heatmap.pdf"), width = 20, height = 7)
for(k in 1:length(p_list)){
    print(p_list[[k]])
}
dev.off()

p_list <- list()
index <- 1
for(k in c(3:12)){
    cat(k, "\n")
    p_list[[index]] <- dittoHeatmap(spe, subregion_genes$D1_islands, group.by = paste0("precast_k", k), 
     order.by = paste0("precast_k", k), scale = "none", heatmap.colors = colorRampPalette(c("white", "red"))(50))
    index <- index + 1
}
pdf(file.path(plot_dir, "D1_islands_marker_expression_heatmap.pdf"), width = 20, height = 7)
for(k in 1:length(p_list)){
    print(p_list[[k]])
}
dev.off()


################################################################################
#   Evaluate if PRECAST cluster assignments were shared for overlapping spots
################################################################################

colData(spe)$overlap_slide <- as.character(colData(spe)$overlap_slide)
colData(spe)$sample_id_original <- as.character(colData(spe)$sample_id_original)
colData(spe)$precast_k2 <- as.factor(colData(spe)$precast_k2)
#   Given a character(1) 'this_key' (a unique barcode in 'col_data$key'), and
#   a tibble of colData 'col_data', return a character(1) giving the key of the
#   spot on a different capture area overlapping the spot given by 'this_key',
#   or return "" if no such spot exists
get_overlapping_key <- function(this_key, col_data) {
    i <- match(this_key, col_data$key)
    x <- col_data |>
        filter(
            sample_id_original == col_data[[i, "overlap_slide"]],
            array_row_transformed == col_data[[i, "array_row_transformed"]],
            array_col_transformed == col_data[[i, "array_col_transformed"]]
        ) |>
        #   Generally, we only had at most 2 capture areas overlapping, but due
        #   to rounding complexities there occasionally can be more than one
        #   spot for the same capture area and array coordinates (hence taking
        #   the first row)
        slice_head(n = 1) |>
        pull(key)

    if (length(x) == 0) {
        return("")
    } else {
        return(x)
    }
}

col_data_full <- colData(spe) |>
    as_tibble()

#   Here, just take the non-excluded spots that overlap another capture area
col_data_small <- col_data_full |>
    filter(!is.na(overlap_slide), overlap_slide != "None", !exclude_overlapping) |>
    select(matches("key|^sample_id$|^precast_k")) |>
    #   Effectively removes spots where PRECAST has not assigned a cluster
    #   identity
    na.omit()

#   Add the key of each overlapping spot, and remove spots where a key doesn't
#   exist
col_data_small <- col_data_small |>
    mutate(
        overlap_key = map_chr(
            col_data_small$key, get_overlapping_key,
            col_data = col_data_full
        )
    ) |>
    filter(overlap_key != "") |>
    #   Append '_original' to the PRECAST clustering results, to signify these
    #   are the results for the non-excluded spots
    rename_with(~ ifelse(grepl("^precast_k", .x), paste0(.x, "_original"), .x))

results_df <- col_data_small |>
    #   Add PRECAST results from the corresponding excluded overlapping spots,
    #   adding '_overlap' to the colnames to differentiate them
    cbind(
        col_data_full[match(col_data_small$overlap_key, col_data_full$key), ] |>
            select(matches("^precast_k")) |>
            rename_with(~ paste0(.x, "_overlap"))
    ) |>
    as_tibble() |>
    #   The colnames now include info about both the capture area (original vs.
    #   overlap) and value of k. Pivot longer to break into 3 columns: capture
    #   area, k, and the cluster identity (assignment)
    pivot_longer(
        cols = matches("^precast_k[0-9]+_"),
        names_to = c("k", "capture_area"),
        names_pattern = "^precast_k([0-9]+)_(overlap|original)$",
        values_to = "cluster_assignment"
    ) |>
    #   Pivot wider so we can compare the cluster identities for the original
    #   and overlapping spots side by side
    pivot_wider(
        names_from = "capture_area", values_from = "cluster_assignment"
    ) |>
    #   Fix a data type
    mutate(k = factor(as.integer(k))) |>
    #   Rarely, an overlapping spot but not the original may have been dropped
    #   as input to PRECAST; simply drop rows with this situation
    na.omit()

#   Plot the proportion of (unique) spots that match their overlapping spot.
#   Each boxplot contains all donors, and we split by value of k
p <- results_df |>
    group_by(sample_id, k) |>
    summarize(match_rate = mean(original == overlap)) |>
    ungroup() |>
    ggplot() +
    geom_boxplot(aes(x = k, y = match_rate), outlier.shape = NA) +
    labs(x = "PRECAST k Value", y = "Proportion of Matches", color = "Donor") +
    geom_jitter(aes(x = k, y = match_rate, color = sample_id)) +
    theme_bw(base_size = 24)

pdf(
    file.path(plot_dir, "cluster_assignment_at_overlaps.pdf"),
    width = 15, height = 7
)
print(p)
dev.off()

#   After some interactive exploration, the proportion of overlapping spots
#   belonging to white matter seems to be the best predictor I found of
#   cluster-assignment mismatches at overlap

#   First, prove the relationship exists at the spot level
message("Among overlapping spots, correlation between matching cluster assignment at")
message("k = 2 and at least one spot in the pair belonging to white matter:")
results_df |>
    filter(k == 2) |>
    summarize(
        wm_vs_match = cor(original == overlap, (original == 2) | (overlap == 2))
    ) |>
    pull(wm_vs_match) |>
    round(2) |>
    message()

#   Next, visually show how the effect manifests at the donor level
p = results_df |>
    filter(k == 2) |>
    group_by(sample_id) |>
    summarize(
        match_rate = mean(original == overlap),
        prop_white = mean((original == 2) | (overlap == 2))
    ) |>
    ungroup() |>
    ggplot() +
        geom_point(aes(x = prop_white, y = match_rate, color = sample_id)) +
        labs(
            x = 'Prop. Pairs w/ at Least 1 WM Spot',
            y = 'Prop. Pairs w/ Matching k = 2 Assignment',
            color = 'Donor'
        ) +
        theme_bw(base_size = 20)
pdf(file.path(plot_dir, "WM_vs_match_rate.pdf"))
print(p)
dev.off()

### For the remainder of the analysis, we can exclude overlapping spots


################################################################################
#   Compare PRECAST clusters with marker gene expression of suspected biological
#   regions (e.g. WM/GM, D1 islands, shell, core)
################################################################################

#   We note that 
b <- Sys.time()
assay_subset <- Matrix(
    assays(spe)$logcounts[match(unlist(subregion_genes), rowData(spe)$gene_name), ],
    sparse = TRUE,
    nrow = length(unlist(subregion_genes)),
    dimnames = list(unlist(subregion_genes), colnames(spe))
)
a <- assay_subset[subregion_genes[[names(subregion_genes)[[1]]]], ]
gene_z <- colMeans((a - rowMeans(a)) / rowSds(a), na.rm = TRUE)
Sys.time() - b

#   Bring in the subset of the logcounts that we'll use, as this seems to be
#   significantly faster than repeatedly performing HDF5-backed, chunked
#   operations
nClust <- 10
cor_mat <- matrix(0, nrow = length(subregion_genes), ncol = nClust)
rownames(cor_mat) <- names(subregion_genes)
for (subregion in names(subregion_genes)) {
    #   For each spot, average expression Z-scores across the subregion.
    #   The idea is that the combined expression of several genes
    #   should form a more precise marker of the subregion than one gene
    a <- assay_subset[subregion_genes[[subregion]], ]
    gene_z <- colMeans((a - rowMeans(a)) / rowSds(a), na.rm = TRUE)

    #   Now compute correlation between cluster identity (in a boolean sense--
    #   belonging or not belonging in the cluster) and the averaged Z-score
    for (cluster_identity in 1:nClust) {
        cor_mat[match(subregion, names(subregion_genes)), cluster_identity] <-
            cor(
                gene_z, spe$precast_k10 == cluster_identity,
                use = "complete.obs"
            )
    }
}

message("Correlation of subregion markers with clusters in k=10:")
print(cor_mat)

p <- pheatmap(cor_mat)
pdf(file.path(plot_dir, "Correlation_of_gene_scores_with_clustering_k10.pdf"))
print(p)
dev.off()

session_info()
