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

/*
 * Requirements for graph class
 * 1. Extract neighbors (out degree)
 * 2. Compute function of neighbor weights and multiply times dense
 * 3. Compute conductivity per-node
 */
template<typename T>
class Graph {
private:
  /* In is major->minor,
   * Out is minor->major 
   * so in col major, row->col is outedge, col->row is inedge
   */
  const Eigen::SparseCompressedBase<T> &m_adj;
  /* In and out weighted degrees (row/col sums) */
  Eigen::ArrayXd m_iwdegree, m_owdegree;
  /* In and out integral degrees (nonzero row/col sums) */
  Eigen::ArrayXi m_iidegree, m_oidegree;
public:
  Graph(const Eigen::SparseCompressedBase<T> &adj,
        Eigen::ArrayXd owdegree, Eigen::ArrayXd iwdegree,
        Eigen::ArrayXi oidegree, Eigen::ArrayXi iidegree) : m_adj(adj),
                                                            m_iwdegree(iwdegree), m_owdegree(owdegree),
                                                            m_iidegree(iidegree), m_oidegree(oidegree) {}
  Graph(const Eigen::SparseCompressedBase<T> &adj) : m_adj(adj),
                                                     m_iwdegree(Eigen::ArrayXd::Zero(T::IsRowMajor ? adj.rows() : adj.cols())),
                                                     m_owdegree(Eigen::ArrayXd::Zero(T::IsRowMajor ? adj.cols() : adj.rows())),
                                                     m_iidegree(Eigen::ArrayXi::Zero(T::IsRowMajor ? adj.rows() : adj.cols())),
                                                     m_oidegree(Eigen::ArrayXi::Zero(T::IsRowMajor ? adj.cols() : adj.rows()))
  { /* Count the in and out degrees (weighted and unweighted) */
#if defined(_OPENMP)
#pragma omp parallel
#endif
  	{
        	Eigen::ArrayXd in_w_degree = Eigen::ArrayXd::Zero(m_iwdegree.size());
                Eigen::ArrayXi in_i_degree = Eigen::ArrayXd::Zero(m_iidegree.size());
#if defined(_OPENMP)
#pragma omp for nowait
#endif
                for (Eigen::Index outerIndex = 0; outerIndex < m_owdegree.size(); outerIndex++) {
                	int local_outer_deg_i = 0;
                        double local_outer_deg_w = 0;
                        for (typename T::InnerIterator it(m_adj, outerIndex); it; ++it) {
                        	Eigen::Index innerIndex = T::IsRowMajor ? it.col() : it.row();
                                if (innerIndex != outerIndex) {
                                	local_outer_deg_i++;
                                        local_outer_deg_w += it.value();
                                        in_i_degree(innerIndex)++;
                                        in_w_degree(innerIndex) += it.value();
                                }
                        }
                        m_oidegree(outerIndex) = local_outer_deg_i;
                        m_owdegree(outerIndex) = local_outer_deg_w;
                }
#if defined(_OPENMP)
#pragma omp critical
#endif
                {
                	m_iidegree += in_i_degree;
                        m_iwdegree += in_w_degree;
                }
        }
  }
  const Eigen::ArrayXd& degree_weighted(bool out=true) const { return out ? m_owdegree : m_iwdegree; }
  const Eigen::ArrayXi& degree_unweighted(bool out=true) const { return out ? m_oidegree : m_iidegree; }
  // std::vector<std::pair<Eigen::Index, T> > neighbors(Eigen::Index outerIndex, bool weighted=true) const {
    
  // }
};

template<typename T>
const Eigen::ArrayXd& degree(const Eigen::SparseCompressedBase<T> &mat)
{
	Graph G(mat);
        return G.degree_weighted();
}

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
