suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("BayesSpace")
    library("RColorBrewer")
    library("ggplot2")
    library("gridExtra")
    library("Polychrome")
})

set.seed(20230712)

# Specify directories
dir_plots <- here("plots", "05_harmony_BayesSpace", "05-BayesSpace_k_search")
spe_in = file.path(here("processed-data", "05_harmony_BayesSpace", "04-preprocess_and_harmony") ,"spe_harmony.rds")
dir_rdata <- here("processed-data", "05_harmony_BayesSpace", "05-BayesSpace_k_search")
spe <- readRDS(spe_in)
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)
colData(spe)$row <- spe$array_row
colData(spe)$col <- spe$array_col

# Read in clustering results
out_path <- here(dir_rdata, "BayesSpace_harmony_k%s", "clusters.csv")
results_list <- lapply(c(2:28), function(K){
     K_nice <- sprintf("%02d", K)
      K_nice <- sprintf("%02d", K)
    df <- sprintf(out_path, K_nice) |>
        read.csv() |>
        as_tibble() |>
        select(c(key, cluster)) |>
        mutate(k = K)
    df
})
   
bayesSpace_results <- do.call(rbind, results_list) |>
    pivot_wider(
        values_from = cluster, names_from = k, names_prefix = "bayesSpace_k"
    )

bayesSpace_results <- data.frame(bayesSpace_results)
for(i in c(2:dim(bayesSpace_results)[2])){
    bayesSpace_results[ ,i] <- factor(as.character(bayesSpace_results[ ,i]), levels = c(1:i))
}

spe$key <- paste0(spe$key, "_", spe$sample_id)
temp <- colnames(spe)
colData(spe) <- colData(spe) |>
    as_tibble() |>
    left_join(bayesSpace_results, by = "key") |>
    DataFrame()
colnames(spe) <- temp

# Check spot plots
for (k in 2:28) {
    plot_list <- list()
    for (donor in levels(spe$sample_id)) {
        plot_list[[donor]] <- spot_plot(
                spe,
                sample_id = donor, var_name = paste0("bayesSpace_k", k),
                is_discrete = TRUE, spatial = TRUE
            ) +
            #   Increase size of colored dots in legend
            guides(fill = guide_legend(override.aes = list(size = 5)))
    }

    pdf(file.path(dir_plots, sprintf("k%s_all_samples.pdf", k)))
    print(plot_list)
    dev.off()
}

