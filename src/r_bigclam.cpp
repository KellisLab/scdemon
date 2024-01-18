
#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include "bigclam.hpp"

// [[Rcpp::export]]
Eigen::MatrixXd r_bigclam_init(const Eigen::Map<Eigen::SparseMatrix<double> > &X,
			  int n_comm)
{
  return bigclam_init(X, n_comm);
}


// [[Rcpp::export]]
Eigen::ArrayXXd r_bigclam(const Eigen::Map<Eigen::SparseMatrix<double> >&C,
			 int n_comm, int max_iter=1000)
{
  return bigclam(C, n_comm, max_iter);
}
