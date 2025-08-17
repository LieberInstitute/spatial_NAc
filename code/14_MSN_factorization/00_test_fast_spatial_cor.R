# Load required packages
rm(list = ls())
library(Matrix)
library(Rcpp)
library(MERINGUE)

# Compile the C++ function
sourceCpp("spatial_fast.cpp")

# Define your wrapper
spatialCrossCorMatrix_fast <- function(mat, weight) {
 # Expect: mat = genes x cells
  #         weight = cells x cells

  if (nrow(weight) != ncol(weight)) {
    stop("'weight' must be a square matrix (cells × cells)")
  }

  if (ncol(mat) != nrow(weight)) {
    stop("'mat' must be genes × cells, where ncol(mat) == nrow(weight)")
  }

  if (!all(rownames(weight) %in% colnames(mat))) {
    stop("Column names of 'mat' (cells) must match rownames of 'weight'")
  }

  # Reorder columns of mat to match weight
  mat <- mat[, rownames(weight)]

  # Convert weight to sparse if necessary
  if (!inherits(weight, "dgCMatrix")) {
    weight <- as(weight, "dgCMatrix")
  }

  # Compute spatial cross-correlation
  scc <- spatialCrossCorMatrix_sparse_C(as.matrix(mat), weight)

  # Name output
  rownames(scc) <- rownames(mat)
  colnames(scc) <- rownames(mat)

  return(scc)
}

# Run test
data(mOB)
pos <- mOB$pos
cd <- mOB$counts

# Remove poor datasets and genes
counts <- cleanCounts(counts = cd, 
                      min.reads = 100, 
                      min.lib.size = 100, 
                      plot=TRUE,
                      verbose=TRUE)

pos <- pos[colnames(counts),]

# CPM normalize
mat <- normalizeCounts(counts = counts, 
                       log=FALSE,
                       verbose=TRUE)

# Sample 2000 genes for demo purposes only to minimize runtime for demo only
set.seed(0)
test <- sample(rownames(mat), 2000)

w <- getSpatialNeighbors(pos, filterDist = 2.5)


start_time <- Sys.time()
I <- getSpatialPatterns(mat[test,], w)
end_time <- Sys.time()
print(end_time - start_time)

results.filter <- filterSpatialPatterns(mat = mat[test,],
                                        I = I,
                                        w = w,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = 0.05,
                                        verbose = TRUE)

scc <- spatialCrossCorMatrix(mat = as.matrix(mat[results.filter,]), 
                             weight = w)

scc_2 <- spatialCrossCorMatrix_fast(mat = as.matrix(mat[results.filter,]), 
                             weight = w)

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

compare_scc_matrices(scc, scc_2)