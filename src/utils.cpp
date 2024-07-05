#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MappedSparseMatrix;

Eigen::MatrixXd createBeta(const Eigen::Map<Eigen::MatrixXd>& G) {
  int n = G.cols();
  Eigen::VectorXd colSum = G.colwise().sum();

  for (int i = 0; i < n; ++i) {
    if (colSum(i) != 0.0) {
      colSum(i) = 1.0 / colSum(i);
    }
  }

  Eigen::MatrixXd betaG = colSum.asDiagonal();
  
  return betaG;
}

// [[Rcpp::export()]]
IntegerVector pred_k(const Eigen::Map<Eigen::VectorXd> & x_mean_square_kK,
            const Eigen::MappedSparseMatrix<double> & query_data,
            const Eigen::Map<Eigen::MatrixXd> & x_mean_kK,
            const Eigen::Map<Eigen::VectorXd> & N_k
)
{
  int K = N_k.size();
  int I = query_data.cols();
  IntegerVector result(I);
  
  for (int i = 0; i < I; i++) {
    int min_index = 0;
    double min_value = x_mean_square_kK(min_index) - (query_data.col(i) + x_mean_kK.col(min_index)).squaredNorm() / (N_k(min_index) + 1);
    
    for (int k = 1; k < K; k++) {
      double current_value = x_mean_square_kK(k) - (query_data.col(i) + x_mean_kK.col(k)).squaredNorm() / (N_k(k) + 1);
      if (current_value < min_value) {
        min_index = k;
        min_value = current_value;
      }
    }
    result[i] = min_index + 1;
  }
  return result;
}