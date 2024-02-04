
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
Eigen::VectorXd r_fw_meat(const Eigen::Map<Eigen::MatrixXd> &res,
			     const Eigen::Map<Eigen::MatrixXd> &U,
			     const Eigen::Map<Eigen::MatrixXd> &B,
			     const Eigen::Map<Eigen::MatrixXd> &BPU)
{
	return fw_meat(res, U, B, BPU);
}

// [[Rcpp::export]]
long double r_fw_bread(const Eigen::Map<Eigen::MatrixXd> &X,
			   const Eigen::Map<Eigen::MatrixXd> &U,
			   const Eigen::Map<Eigen::MatrixXd> &B,
			   const Eigen::Map<Eigen::MatrixXd> &BPU)
{
	return fw_bread(X, U, B, BPU);
}

// [[Rcpp::export]]
Eigen::MatrixXd r_cbind(const Eigen::Map<Eigen::MatrixXd> &X,
			const Eigen::Map<Eigen::MatrixXd> &Y)
{
	return cbind(X, Y);
}

// [[Rcpp::export]]
Eigen::ArrayXd r_robust_se_X(const Eigen::Index &x_idx,
			     const Eigen::Map<Eigen::MatrixXd> &Y)
{
	return robust_se_X(x_idx, Y, 1e-300);
}

// [[Rcpp::export]]
Eigen::ArrayXd r_robust_se_Xfull(const Eigen::Index &x_idx,
                                 const Eigen::Map<Eigen::MatrixXd> &U,
                                 const Eigen::Map<Eigen::MatrixXd> &V,
                                 const Eigen::Index block_size)
{
	return robust_se_Xfull(x_idx, U, V, 1e-300, block_size);
}
// [[Rcpp::export]]
Eigen::SparseMatrix<double> r_robust_se(const Eigen::Map<Eigen::MatrixXd> &Y,
					double t_cutoff=2, bool abs_cutoff=false)
{
	return robust_se(Y, NULL, NULL, 1e-300, t_cutoff, abs_cutoff);
}
			   


// [[Rcpp::export]]
Eigen::SparseMatrix<double> r_robust_se_p(const Eigen::Map<Eigen::MatrixXd> &Y,
					  const Eigen::Map<Eigen::ArrayXd> &dof,
					  double nominal_p_cutoff=0.05, bool abs_cutoff=false)
{
  #if defined(_USE_GSL)
	return robust_se_pvalue(Y, dof, NULL, NULL, nominal_p_cutoff, abs_cutoff, 1e-300);
  #else
	throw Rcpp::exception("Unimplemented due to lack of GSL");
  #endif
}
