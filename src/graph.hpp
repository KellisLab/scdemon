#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <vector>
#include <unordered_set>
#include <utility>
#include <algorithm>

#include "implprogress.hpp"
template<typename T>
std::vector<std::pair<Eigen::Index, double> > graph_neighbors(const Eigen::SparseCompressedBase<T> &mat,
							      Eigen::Index outerIndex,
							      bool weighted=true)
{
        std::vector<std::pair<Eigen::Index, double> > out;
	for (typename T::InnerIterator it(mat, outerIndex); it; ++it) {
		out.emplace_back(T::IsRowMajor ? it.col() : it.row(),
				 weighted ? 1. / it.value() : 1.);
	}
	return out;
}

/* Since Eigen 3.4 is not supported in Rcpp, write own iterator */
template<typename T>
Eigen::MatrixXd matrix_slice_rows(const Eigen::DenseBase<T> &mat,
				  const std::vector<std::pair<Eigen::Index, double> > &indices)
{
        Eigen::MatrixXd out(indices.size(), mat.cols());
	for (std::pair<Eigen::Index, double> i : indices) {
	  	Eigen::ArrayXd row = mat.row(i.first).array() * i.second;
	        out << row.matrix();
	}
	return out;
}

template<typename T>
double graph_conductance_row(const Eigen::SparseCompressedBase<T> &mat, Eigen::Index i, const Eigen::ArrayXd& degree, double deg_sum, int min_deg=5)
{
	std::unordered_set<Eigen::Index> subgraph = {i};
	for (typename T::InnerIterator it(mat, i); it; ++it) {
		Eigen::Index j = T::IsRowMajor ? it.col() : it.row();
		if (i != j) { 
			subgraph.insert(j);
		}
	}
	if (subgraph.size() < min_deg) {
		return 1.;
	}
	double subgraph_sum = 0;
	double subgraph_degree_sum = 0;
	for (Eigen::Index j : subgraph) {
		// For each neighbor, add the degrees
		subgraph_degree_sum += degree(j);
		for (typename T::InnerIterator it(mat, j); it; ++it) {
			Eigen::Index k = T::IsRowMajor ? it.col() : it.row();
			if (j == k) {
				continue;
			} else if (subgraph.find(k) != subgraph.end()) {
				subgraph_sum += it.value();
			}
		}
	}
	/*
	 * Numerator is the edge crossings:
	 * Denom is deg sum and its complement
	 */
	double numer = subgraph_degree_sum - subgraph_sum;
	double denom = std::min(subgraph_degree_sum, deg_sum - subgraph_degree_sum);
	return numer / denom;
}

template<typename T>
Eigen::ArrayXd graph_degree(const Eigen::SparseCompressedBase<T> &mat, bool diag=false)
{
	Eigen::Index size;
	if (T::IsRowMajor) {
		size = mat.rows();
	} else {
		size = mat.cols();
	}
	Eigen::ArrayXd degree = Eigen::ArrayXd::Zero(size);
#if defined(_OPENMP)
#pragma omp parallel for
#endif
	for (Eigen::Index i = 0; i < size; i++) {
		double cur = 0;
		for (typename T::InnerIterator it(mat, i); it; ++it) {
			Eigen::Index j = T::IsRowMajor ? it.col() : it.row();
			if (diag || (i != j)) {
				cur += it.value();
			}
			degree(i) = cur;
		}
	}
	return degree;
}
template<typename T>
Eigen::ArrayXd graph_conductance(const Eigen::SparseCompressedBase<T> &mat, int min_deg=5)
{
	Eigen::ArrayXd degree = graph_degree(mat, false);
	const double deg_sum = degree.sum();
	ImplProgress p(mat.cols());
	Eigen::ArrayXd conductance = Eigen::ArrayXd::Ones(degree.size());
#if defined(_OPENMP)
#pragma omp parallel for
#endif
	for (Eigen::Index i = 0; i < mat.cols(); i++) {
	/*
	 * For each index, find the subgraph defined
	 * by the immediate neighbors of row i
	 */
		if (!p.check_abort()) {
			p.increment();
			double cond = graph_conductance_row(mat, i, degree, deg_sum, min_deg);
			conductance(i) = cond;
		}
	}
	return conductance;
}
#endif
