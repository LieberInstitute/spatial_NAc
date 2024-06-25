library(here)
library(PRECAST)
library(HDF5Array)
library(Seurat)
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

nRandom_starts <- 4
K <- 3:15
resDir <- here("processed-data", "10_precast", "nnSVG_precast")

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
    p_early_iterations [[i]] <- ggplot(convergence_df[convergence_df$iteration > 10 & convergence_df$iteration < 50, ], aes(x = iteration, y = loglik, color = random_start)) +
    geom_point() + geom_line() + theme_classic() + ggtitle(paste0("K=", convergence_df$nClusters[1]))
}
plot_grid(p_early_iterations[[1]], p_early_iterations[[2]], p_early_iterations[[3]], p_early_iterations[[4]], 
p_early_iterations[[5]], p_early_iterations[[6]], nrow = 3)

plot_grid(p_early_iterations[[7]], p_early_iterations[[8]],
p_early_iterations[[9]], p_early_iterations[[10]], p_early_iterations[[11]], p_early_iterations[[12]], ncol = 2)

p_late_iterations <- list()
for(i in c(1:length(convergence_results))){
    convergence_df <- convergence_results[[i]]
    p_late_iterations [[i]] <- ggplot(convergence_df[convergence_df$iteration > 100, ], aes(x = iteration, y = loglik, color = random_start)) +
    geom_point() + geom_line() + theme_classic() + ggtitle(paste0("K=", convergence_df$nClusters[1]))
}

plot_grid(p_late_iterations[[1]], p_late_iterations[[2]], p_late_iterations[[3]], p_late_iterations[[4]], 
p_late_iterations[[5]], p_late_iterations[[6]], nrow = 3)

plot_grid(p_late_iterations[[7]], p_late_iterations[[8]], p_late_iterations[[9]], p_late_iterations[[10]], 
p_late_iterations[[11]], p_late_iterations[[12]], nrow = 3)

summary_list <- list()
for(i in c(1:length(convergence_results))){
    convergence_df <- convergence_results[[i]]
    max_ll_1 <- max(convergence_df[convergence_df$random_start == "1", "loglik"])
    max_ll_2 <-max(convergence_df[convergence_df$random_start == "2", "loglik"])
    max_ll_3 <-max(convergence_df[convergence_df$random_start == "3", "loglik"])
    summary_df <- data.frame("random_start" = c(1, 2, 3), 
    "max_ll" = c(max_ll_1, max_ll_2, max_ll_3))
    summary_df$cluster <- K[i]
    summary_df$random_start <- factor(summary_df$random_start, levels = c(1:nRandom_starts))
    summary_list[[i]] <- summary_df 
}

p_final_ll <- list()
for(i in c(1:length(summary_list))){
    p_final_ll[[i]] <- ggplot(summary_list[[i]], aes(x = random_start, y = max_ll, color = random_start)) + geom_point() + theme_bw() +
ggtitle(paste0("K=", summary_list[[i]]$cluster[1])) + xlab("Random start") + ylab("Final log-lik")
}

plot_grid(p_final_ll[[1]], p_final_ll[[2]], p_final_ll[[3]], p_final_ll[[4]], 
p_final_ll[[5]], p_final_ll[[6]], nrow = 3)

plot_grid(p_final_ll[[7]], p_final_ll[[8]], p_final_ll[[9]], p_final_ll[[10]], 
p_final_ll[[11]], p_final_ll[[12]], nrow = 3)