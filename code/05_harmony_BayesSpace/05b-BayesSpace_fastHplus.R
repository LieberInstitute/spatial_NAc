library(fasthplus)
library(SpatialExperiment)
library(here)
library(spatialLIBD)

## Choose k
k <- as.numeric(
    #   Only one of these environment variables will be defined, so grab the
    #   defined one (handle SGE or SLURM)
    paste0(Sys.getenv("SLURM_ARRAY_TASK_ID"), Sys.getenv("SGE_TASK_ID"))
)
k_nice <- sprintf("%02d", k)

# Read spe
spe_in = file.path(here("processed-data", "05_harmony_BayesSpace", "04-preprocess_and_harmony") ,"spe_harmony.rds")
spe <- readRDS(spe_in)

# Read in BayesSpace results
dir_rdata <- here("processed-data", "05_harmony_BayesSpace", "05-BayesSpace_k_search")
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

## remove white matter
dim(spe)
# [1]  28916 113927
spe <- spe[, -which(colData(spe)$bayesSpace_harmony_02 == 1)]
dim(spe)
