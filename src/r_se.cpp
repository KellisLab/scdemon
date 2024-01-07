
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
Eigen::VectorXd r_adjust_res(const Eigen::Map<Eigen::MatrixXd> &res,
			     const Eigen::Map<Eigen::MatrixXd> &U,
			     const Eigen::Map<Eigen::MatrixXd> &B,
			     const Eigen::Map<Eigen::MatrixXd> &BPU)
{
        return adjust_res(res, U, B, BPU);
}

// [[Rcpp::export]]
Eigen::VectorXd r_hc0_se_Xvec(const Eigen::Map<Eigen::MatrixXd> &X,
			      const Eigen::Map<Eigen::MatrixXd> &Vs,
			      const Eigen::Map<Eigen::MatrixXd> &U,
			      const Eigen::Map<Eigen::MatrixXd> &TUU,
			      const Eigen::Map<Eigen::MatrixXd> &B,
			      const Eigen::Map<Eigen::MatrixXd> &BPU) {
        return hc0_se_Xvec(X, Vs, U, TUU, B, BPU);
}
