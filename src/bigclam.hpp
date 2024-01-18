#ifndef BIGCLAM_HPP
#define BIGCLAM_HPP
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <algorithm>
#include <vector>

#include "graph.hpp"
template<typename T>
Eigen::MatrixXd bigclam_init(const Eigen::SparseMatrixBase<T> &C, int n_comm)
{
        Eigen::ArrayXXd F = Eigen::ArrayXXd::Random(std::min(C.rows(), C.cols()),
						    n_comm);
	return ((F + 1.) * 0.5).matrix();
}

template<typename TV, typename TN, typename TF>
Eigen::ArrayXd bigclam_gradient(const Eigen::MatrixBase<TV> &node,
				const Eigen::MatrixBase<TN> &neighbors,
				const Eigen::ArrayBase<TF> &Fsum0)
{
  Eigen::ArrayXd raw = (node * neighbors.transpose()).array().min(15).max(-15).eval();
  Eigen::ArrayXd scores = (-raw).exp() / (1 - (-raw).exp());
  Eigen::ArrayXd neighbor_grad = (scores.matrix() * neighbors).array().eval();
  Eigen::ArrayXd no_neigh_grad = (Fsum0 - node.array() - neighbors.array().colwise().sum()).eval();
  return neighbor_grad - no_neigh_grad;
}

// bg edge: \epsilon = 2|E|/(|V||V-1|)
// todo find out why F increases with # of iterations
#include <iostream>
template<typename T>
Eigen::ArrayXXd bigclam(const Eigen::SparseCompressedBase<T> &C,
			int n_comm, int max_iter=1000,
			double conv_threshold=1e-3, double eta=0.005,
			double min_value=0.00001, double max_value=20)
{
        Eigen::MatrixXd F = bigclam_init(C, n_comm);
	Eigen::ArrayXd Fsum0 = F.colwise().sum();
	std::cout << "Fsum0 size: " <<  Fsum0.size() << " data: " << Fsum0 << std::endl;
	for (int iter = 0; iter < max_iter; iter++) {
	  std::cout << "Iter: " << iter << std::endl;
	  for (Eigen::Index i = 0; i < C.cols(); i++) {
	    std::vector<Eigen::Index> neb_idx = graph_neighbors(C, i);
	    //std::cout << i << " Neighbor: " << neb_idx.size() << std::endl;
	    if (neb_idx.size() > 0) {
	      Eigen::MatrixXd F_neighbors = matrix_slice_rows(F, neb_idx);
	      Eigen::ArrayXd F_node = F.row(i);
	      Eigen::ArrayXd gradient = bigclam_gradient(F_node.matrix(), F_neighbors, Fsum0);
	      F.row(i) = F.row(i).array() + eta * gradient;
	      F.row(i) = F.row(i).array().max(min_value).min(max_value);
	      Fsum0 = Fsum0 - F_node + F.row(i).array();
	    }
	  }
	}
	return F.array();
}
#endif
