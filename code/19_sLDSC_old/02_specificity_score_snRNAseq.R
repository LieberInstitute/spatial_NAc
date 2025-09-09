# Load libraries
library(here)

is_top_10_percent <- function(column) {
    top_10_percent_threshold <- quantile(column, 0.9)
    as.numeric(column >= top_10_percent_threshold)
}

dat <- read.table(here::here("processed-data", "19_sLDSC_old","snRNA_seq","input_files", "snRNA_aggregated_de.tsv"),header=T)

res <- apply(dat, 2, is_top_10_percent)
rownames(res) <- rownames(dat)
write.csv(res,here::here("processed-data", "19_sLDSC_old", "snRNA_seq", "input_files", "snRNA_score.csv"))