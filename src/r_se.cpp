
#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include "se.hpp"

// [[Rcpp::export]]
Eigen::SparseMatrix<double> r_robust_se(const Eigen::Map<Eigen::MatrixXd> &U,
					const Eigen::Map<Eigen::MatrixXd> &V,
					double lambda=1e-10,
					double t_cutoff=2, bool abs_cutoff=false)
{
	return intra_robust_se(U, V, lambda, t_cutoff, abs_cutoff);
}
