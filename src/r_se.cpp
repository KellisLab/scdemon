
#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include "se.hpp"

  
// [[Rcpp::export]]
Eigen::MatrixXd r_ols_beta(const Eigen::Map<Eigen::MatrixXd> &X,
			   const Eigen::Map<Eigen::MatrixXd> &Y)
{
	return ols_beta(X, Y);
}

// [[Rcpp::export]]
Eigen::MatrixXd r_ols_resid(const Eigen::Map<Eigen::MatrixXd> &X,
			    const Eigen::Map<Eigen::MatrixXd> &Y,
			    const Eigen::Map<Eigen::MatrixXd> &beta)
{
	return ols_resid(X, Y, beta);
}

// [[Rcpp::export]]
Eigen::ArrayXd r_robust_se_X(const Eigen::Map<Eigen::MatrixXd> &X,
			     const Eigen::Map<Eigen::MatrixXd> &Y)
{
	return robust_se_X(X, Y, 1e-300);
}


// [[Rcpp::export]]
Eigen::SparseMatrix<double> r_robust_se(const Eigen::Map<Eigen::MatrixXd> &X,
					const Eigen::Map<Eigen::MatrixXd> &Y,
					double t_cutoff=2, bool abs_cutoff=false)
{
	return robust_se(X, Y, 1e-300, t_cutoff, abs_cutoff);
}
