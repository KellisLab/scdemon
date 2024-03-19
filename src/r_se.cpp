
#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include "se.hpp"

  
// [[Rcpp::export]]
Eigen::MatrixXd r_ols_beta(const Eigen::Map<Eigen::MatrixXd> &X,
			   const Eigen::Map<Eigen::MatrixXd> &Y,
                           double lambda)
{
	return ols_beta(X, Y, lambda);
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
			     const Eigen::Map<Eigen::MatrixXd> &Y,
			     double lambda=1e-10)
{
	return robust_se_X(X, Y, lambda, 1e-300);
}

// [[Rcpp::export]]
Eigen::ArrayXd r_robust_se_Y(const Eigen::Map<Eigen::MatrixXd> &X,
			     const Eigen::Map<Eigen::MatrixXd> &Y,
			     double lambda=1e-10)
{
	return robust_se_Y(X, Y, lambda, 1e-300);
}

// [[Rcpp::export]]
Eigen::MatrixXd r_ols_beta_L(const Eigen::Map<Eigen::VectorXd> &X,
                             const Eigen::Map<Eigen::MatrixXd> &Y,
                             const Eigen::Map<Eigen::ArrayXd> &lambda)
{
	return ols_beta_L(X.eval(), X.squaredNorm(), Y, lambda);
}
// [[Rcpp::export]]
Eigen::ArrayXd r_robust_se_L(const Eigen::Map<Eigen::MatrixXd> &X,
			     const Eigen::Map<Eigen::MatrixXd> &Y,
                             const Eigen::Map<Eigen::ArrayXd> &lambda)
{
	return robust_se_L(X, Y, lambda, 1e-300);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> r_robust_se(const Eigen::Map<Eigen::MatrixXd> &X,
					const Eigen::Map<Eigen::MatrixXd> &Y,
					double lambda=1e-10,
					double t_cutoff=2, bool abs_cutoff=false)
{
	return robust_se(X, Y, lambda, 1e-300, t_cutoff, abs_cutoff);
}
