rm(list = ls())
library(Matrix)
library(Rcpp)
library(MERINGUE)

compare_scc_matrices <- function(scc1, scc2, tolerance = 0.05, plot = FALSE) {
  stopifnot(all(dim(scc1) == dim(scc2)))
  stopifnot(all(rownames(scc1) == rownames(scc2)))
  stopifnot(all(colnames(scc1) == colnames(scc2)))

  # Flatten upper triangle (excluding diagonal)
  get_upper <- function(mat) {
    mat[upper.tri(mat)]
  }

  vec1 <- get_upper(scc1)
  vec2 <- get_upper(scc2)

  # Sign agreement
  sign_match <- sign(vec1) == sign(vec2)
  sign_agreement_rate <- mean(sign_match)

  # Magnitude difference
  abs_diff <- abs(vec1 - vec2)
  rel_diff <- abs_diff / (abs(vec1) + 1e-10)  # avoid division by zero

  # Ordering agreement: Spearman rank
  spearman_cor <- cor(vec1, vec2, method = "spearman")

  # Summary
  cat("SCC Matrix Comparison:\n")
  cat("Same dimensions:", dim(scc1)[1], "genes\n")
  cat("Sign agreement rate:", round(sign_agreement_rate * 100, 2), "%\n")
  cat("Mean absolute difference:", signif(mean(abs_diff), 4), "\n")
  cat("Max absolute difference:", signif(max(abs_diff), 4), "\n")
  cat("Mean relative difference:", signif(mean(rel_diff), 4), "\n")
  cat("Spearman rank correlation:", signif(spearman_cor, 4), "\n")

  if (plot) {
    op <- par(mfrow = c(1, 2))
    smoothScatter(vec1, vec2,
                  xlab = "Original SCC",
                  ylab = "Optimized SCC",
                  main = "SCC Value Agreement")
    abline(0, 1, col = "red", lty = 2)

    hist(rel_diff, breaks = 50,
         main = "Relative Differences",
         xlab = "Relative difference")
    par(op)
  }

  invisible(list(
    sign_agreement = sign_agreement_rate,
    spearman = spearman_cor,
    abs_diff = abs_diff,
    rel_diff = rel_diff
  ))
}

sample_ids <- c("Br2743", "Br3942", "Br6423", "Br6432", "Br6471", "Br8492", "Br8667")
res_dir_old <- here("processed-data", "14_MSN_factorization", "01_meringue", "original_script_scc")
res_dir_new <- here("processed-data", "14_MSN_factorization", "01_meringue")

for(i in c(1:length(sample_ids))){
    print(sample_ids[i])
    scc_old <- readRDS(file.path(res_dir_old, paste0("scc_", sample_ids[i], ".rds")))
    scc_new <- readRDS(file.path(res_dir_new, paste0("scc_", sample_ids[i], ".rds")))
    compare_scc_matrices(scc_old, scc_new)
}
