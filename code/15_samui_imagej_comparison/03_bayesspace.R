library(getopt)
library(sessioninfo)
library(here)
library(SpatialExperiment)
library(BayesSpace)
library(tidyverse)
library(HDF5Array)

# Import command-line parameters
spec <- matrix(
    c(
        c("final_step", "k"),
        c("f", "k"),
        rep("1", 2),
        c("character", "numeric"),
        c(
            "Using array coordinates after either 'imagej' or 'samui",
            "Number of clusters"
        )
    ),
    ncol = 5
)
opt <- getopt(spec)

print("Using the following parameters:")
print(opt)

spe_dir = here('processed-data', '15_samui_imagej_comparison', 'spe')
out_path = here(
    'processed-data', '15_samui_imagej_comparison', 'bayesspace_out',
    sprintf('k%s_%s.csv', opt$k, opt$final_step)
)

set.seed(1)
dir.create(dirname(out_path), showWarnings = FALSE)

spe = loadHDF5SummarizedExperiment(spe_dir)

#   Since BayesSpace considers all donors at once, we need to separate donors
#   spatially. Do this by adding a donor-specific offset to the array row. Note
#   here that 'array_row' is not constrained to have max value 77; we instead
#   find the largest 'array_row' value of any donor, and use it to ensure donors
#   are at least 5 rows apart
offset_row <- as.numeric(factor(spe$sample_id)) *
    (max(spe[[paste0('array_row_', opt$final_step)]]) + 5)

#   Set array coordinates based on final step (ImageJ or Samui)
spe$row = spe[[paste0('array_row_', opt$final_step)]] + offset_row
spe$col = spe[[paste0('array_col_', opt$final_step)]]

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

#   Run main clustering step
message("Running spatialCluster()")
Sys.time()
spe <- spatialCluster(spe, use.dimred = "HARMONY", q = opt$k, nrep = 10000)
Sys.time()

#   Write cluster assigments (along with spot key) to CSV
colData(spe) |>
    as_tibble() |>
    select(key, spatial.cluster) |>
    rename(cluster_assignment = spatial.cluster) |>
    write_csv(out_path)

session_info()

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/
