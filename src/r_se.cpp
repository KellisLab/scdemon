
#include <RcppEigen.h>
#include "se.hpp"

// [[Rcpp::export]]
Eigen::MatrixXd r_ols(const Eigen::VectorXd &X, const Eigen::MatrixXd &Y) {
    return ols(X, Y);
}
