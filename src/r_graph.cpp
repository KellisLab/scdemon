
#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include "graph.hpp"

// [[Rcpp::export]]
Eigen::ArrayXd r_graph_conductance(const Eigen::Map<Eigen::SparseMatrix<double> > &X,
				   int min_deg)
{
	return graph_conductance(X, min_deg);
}

// [[Rcpp::export]]
double r_graph_conductance_row(const Eigen::Map<Eigen::SparseMatrix<double> > &X,
				       Eigen::Index i,
				       const Eigen::ArrayXd& degree,
				       double deg_sum,
				       int min_deg)
{
	return graph_conductance_row(X, i, degree, deg_sum, min_deg);
}
// [[Rcpp::export]]
Eigen::ArrayXd r_graph_degree(const Eigen::Map<Eigen::SparseMatrix<double> > &X,
			      bool diag=false)
{
  return graph_degree(X, diag);
}
