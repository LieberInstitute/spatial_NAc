library(here)
library(PRECAST)
library(HDF5Array)
library(Seurat)
library(mclust)
library(ggplot2)
library(cowplot)
library(dplyr)
library(sessioninfo)
library(tidyverse)
library(Matrix)
library(SpatialExperiment)
library(ggsci)
library(ggpubr)
library(getopt)
library(scran)
library(scater)
library(spatialNAcUtils)
library(spatialLIBD)
library(ggpubr)

nRandom_starts <- 5
K <- 3:15
resDir <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast")
plotDir <- here("plots", "07_spatial_domains", "01_precast", "nnSVG_precast", "examine_convergence")
convergence_results <- lapply(K, function(k){
     cat(k, "\n")
     convergence_list <- lapply(c(1:nRandom_starts), function(rs){
        pre_obj <- readRDS(paste0(resDir, "/random_start_", rs, "/pre_obj_k", k, ".rds"))
        ll_seq <- pre_obj@resList[[1]]$loglik_seq
        iter <- c(1:length(ll_seq))
        convergence_df <- data.frame("iteration" = iter,  "loglik" = ll_seq, 
                                     "random_start" = rs)
        convergence_df <- convergence_df[!convergence_df$loglik == 0, ]
        convergence_df
    })
    convergence_df <- do.call(rbind, convergence_list)
    convergence_df$nClusters <- k
    convergence_df
})

convergence_results <- lapply(convergence_results, function(convergence_df){
    convergence_df$random_start <- factor(convergence_df$random_start, levels = c(1:nRandom_starts))
    convergence_df
})
p_early_iterations <- list()
for(i in c(1:length(convergence_results))){
    convergence_df <- convergence_results[[i]]
    p_early_iterations [[i]] <- ggplot(convergence_df[convergence_df$iteration > 1 & convergence_df$iteration <= 20, ], aes(x = iteration, y = loglik, color = random_start)) +
    geom_point() + geom_line() + theme_classic() + ggtitle(paste0("K=", convergence_df$nClusters[1])) + guides(fill = guide_legend(title = "Random start"), color = guide_legend(title = "Random start"))
}

pdf(file.path(plotDir, "early_iterations.pdf"))
plot_grid(p_early_iterations[[1]], p_early_iterations[[2]], p_early_iterations[[3]], p_early_iterations[[4]], 
p_early_iterations[[5]], p_early_iterations[[6]], nrow = 3)
plot_grid(p_early_iterations[[7]], p_early_iterations[[8]],
p_early_iterations[[9]], p_early_iterations[[10]], p_early_iterations[[11]], p_early_iterations[[12]], ncol = 2)
plot_grid(p_early_iterations[[13]], ncol = 2, nrow = 3)
dev.off()

p_late_iterations <- list()
for(i in c(1:length(convergence_results))){
    convergence_df <- convergence_results[[i]]
    p_late_iterations [[i]] <- ggplot(convergence_df[convergence_df$iteration > 100, ], aes(x = iteration, y = loglik, color = random_start)) +
    geom_point() + geom_line() + theme_classic() + ggtitle(paste0("K=", convergence_df$nClusters[1])) +  guides(fill = guide_legend(title = "Random start"), color = guide_legend(title = "Random start"))
}

pdf(file.path(plotDir, "late_iterations.pdf"))
plot_grid(p_late_iterations[[1]], p_late_iterations[[2]], p_late_iterations[[3]], p_late_iterations[[4]], 
p_late_iterations[[5]], p_late_iterations[[6]], nrow = 3, ncol = 2)
plot_grid(p_late_iterations[[7]], p_late_iterations[[8]], p_late_iterations[[9]], p_late_iterations[[10]], 
p_late_iterations[[11]], p_late_iterations[[12]], nrow = 3, ncol = 2)
plot_grid(p_late_iterations[[13]], nrow = 3, ncol = 2)
dev.off()

summary_list <- list()
for(i in c(1:length(convergence_results))){
    convergence_df <- convergence_results[[i]]
    max_ll_1 <- max(convergence_df[convergence_df$random_start == "1", "loglik"])
    max_ll_2 <-max(convergence_df[convergence_df$random_start == "2", "loglik"])
    max_ll_3 <-max(convergence_df[convergence_df$random_start == "3", "loglik"])
    max_ll_4 <-max(convergence_df[convergence_df$random_start == "4", "loglik"])
    max_ll_5 <-max(convergence_df[convergence_df$random_start == "5", "loglik"])

    summary_df <- data.frame("random_start" = c(1, 2, 3, 4, 5), 
    "max_ll" = c(max_ll_1, max_ll_2, max_ll_3, max_ll_4, max_ll_5))
    summary_df$cluster <- K[i]
    summary_df$random_start <- factor(summary_df$random_start, levels = c(1:nRandom_starts))
    summary_list[[i]] <- summary_df 
}

summary_df <- do.call(rbind, summary_list)

pdf(file.path(plotDir, "loglik_distributions_at_convergence.pdf"), height = 4.5, width = 9)
ggplot(summary_df, aes(x = cluster, y = max_ll, color = random_start, fill = random_start)) + geom_line() + geom_point()  + theme_pubr() +
xlab("Cluster") + ylab("Maximum log-lik") + guides(fill = guide_legend(title = "Random start"), color = guide_legend(title = "Random start"))
dev.off()

spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5"
)
spe <- loadHDF5SummarizedExperiment(spe_dir)

lapply(K, function(k){
    cat(k, "\n")
    cluster_list <- lapply(c(1:nRandom_starts), function(rs){
        precast_results <- read.csv(paste0(resDir, "/random_start_", rs, "/", sprintf("PRECAST_k%s.csv", k)))
        precast_results$nCount <- as.integer(precast_results$nCount)
        precast_results$batch <- factor(precast_results$batch)
        precast_results$cluster <- factor(precast_results$cluster)
        precast_results
})
    stopifnot(all.equal(cluster_list[[1]]$key, cluster_list[[2]]$key))
    stopifnot(all.equal(cluster_list[[2]]$key, cluster_list[[3]]$key))
    stopifnot(all.equal(cluster_list[[3]]$key, cluster_list[[4]]$key))
    stopifnot(all.equal(cluster_list[[4]]$key, cluster_list[[5]]$key))
    
    spe_sub <- spe[ ,colnames(spe) %in% cluster_list[[1]]$key]
    spe_sub <- spe_sub[ ,match(cluster_list[[1]]$key, colnames(spe_sub))]
    spe_sub[[paste0("precast_k", k, "_rs", 1)]] <- cluster_list[[1]]$cluster
    spe_sub[[paste0("precast_k", k, "_rs", 2)]] <- cluster_list[[2]]$cluster
    spe_sub[[paste0("precast_k", k, "_rs", 3)]] <- cluster_list[[3]]$cluster
    spe_sub[[paste0("precast_k", k, "_rs", 4)]] <- cluster_list[[4]]$cluster
    spe_sub[[paste0("precast_k", k, "_rs", 5)]] <- cluster_list[[5]]$cluster

    samples <- as.character(unique(spe_sub$sample_id))
    plot_clust_dir <- paste0(plotDir, "/precast_k", k)
    dir.create(plot_clust_dir, showWarnings = FALSE)
    plot_list <- list()
    for(rs in c(1:nRandom_starts)){
        plot_list[[rs]] <- spot_plot(
            spe_sub,
            sample_id = "Br2720", var_name = paste0("precast_k", k, "_rs", rs),
            is_discrete = TRUE, spatial = TRUE
        ) +
        #   Increase size of colored dots in legend
        guides(fill = guide_legend(override.aes = list(size = 5))) + ggtitle(paste("Random start:", rs))
        }
    pdf(paste0(plot_clust_dir, "/", "Br2720.pdf"), width = 15, height = 12)
    plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], nrow = 2)
    dev.off()

    plot_list <- list()
    for(rs in c(1:nRandom_starts)){
        plot_list[[rs]] <- spot_plot(
            spe_sub,
            sample_id = "Br6522", var_name = paste0("precast_k", k, "_rs", rs),
            is_discrete = TRUE, spatial = TRUE
        ) +
        #   Increase size of colored dots in legend
        guides(fill = guide_legend(override.aes = list(size = 5))) + ggtitle(paste("Random start:", rs))
        }
    pdf(paste0(plot_clust_dir, "/", "Br6522.pdf"), width = 15, height = 10)
    plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], nrow = 2)
    dev.off()
    
    plot_list <- list()
    for(rs in c(1:nRandom_starts)){
        plot_list[[rs]] <- spot_plot(
            spe_sub,
            sample_id = "Br3942", var_name = paste0("precast_k", k, "_rs", rs),
            is_discrete = TRUE, spatial = TRUE
        ) +
        #   Increase size of colored dots in legend
        guides(fill = guide_legend(override.aes = list(size = 5))) + ggtitle(paste("Random start:", rs))
        }
    pdf(paste0(plot_clust_dir, "/", "Br3942.pdf"), width = 15, height = 10)
    plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], nrow = 2)
    dev.off()

    concordance_mat <- matrix(NA, nrow = nRandom_starts, ncol = nRandom_starts)
    rownames(concordance_mat) <- paste0("RS_", c(1:5))
    colnames(concordance_mat) <- paste0("RS_", c(1:5))
    for(i in c(1:nRandom_starts)){
        for(j in c(1:nRandom_starts)){
            if(i != j){
                clust_a <- as.character(cluster_list[[i]]$cluster)
                clust_b <- as.character(cluster_list[[j]]$cluster)
                concordance_mat[i, j] <- adjustedRandIndex(clust_a, clust_b)
            }
        }
    }
    concordance_df <- reshape2::melt(concordance_mat)

    pdf(paste0(plot_clust_dir, "/", "rand_index.pdf"), height = 6, width = 6)
    ggplot(concordance_df, aes(x = Var1, y = Var2, fill = value)) + geom_tile(color = "black") +
    geom_text(aes(label = round(value, 2)), color = "black", size = 4) + scale_fill_gradient(low = "white", high = "red") +
    coord_fixed() + xlab("") + ylab("") + theme_pubr() + theme(axis.ticks = element_blank(), axis.line = element_blank()) 
    dev.off()
})
