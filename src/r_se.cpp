
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
Eigen::VectorXd r_hc0_se_Xvec(const Eigen::Map<Eigen::MatrixXd> &X,
			      const Eigen::Map<Eigen::MatrixXd> &Vs,
			      const Eigen::Map<Eigen::MatrixXd> &U,
			      const Eigen::Map<Eigen::MatrixXd> &B,
			      const Eigen::Map<Eigen::MatrixXd> &BPU) {
        return hc0_se_Xvec(X, Vs, U, B, BPU);
}
