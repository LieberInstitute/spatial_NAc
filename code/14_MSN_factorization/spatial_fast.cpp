// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat spatialCrossCorMatrix_sparse_C(arma::mat sigMat, arma::sp_mat weight) {
  int N = sigMat.n_rows; // genes
  int M = sigMat.n_cols; // cells

  // Step 1: Normalize sparse weight matrix by row
  arma::vec rowSums = arma::vec(arma::sum(weight, 1));
  for (arma::sp_mat::iterator it = weight.begin(); it != weight.end(); ++it) {
    double denom = rowSums(it.row());
    if (denom == 0) denom = 1.0;
    *it = *it / denom;
  }

  // Step 2: Compute W after normalization
  double W = accu(weight);

  // Step 3: Center gene expression (genes x cells)
  arma::mat centered = sigMat.each_col() - arma::mean(sigMat, 1);

  // Step 4: Spatially weighted expression
  arma::mat weighted = centered * weight.t();  // genes × cells

  // Step 5: Cross-correlation
  arma::mat cross = centered * weighted.t();   // genes × genes

  // Step 6: Normalize by variance
  arma::vec variances = sum(square(centered), 1);
  arma::mat denom = sqrt(variances * variances.t());
  arma::mat scc = (M / W) * (cross / denom);

  return scc;
}
