#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <vector>
template<typename T>
std::vector<Eigen::Index> graph_neighbors(const Eigen::SparseCompressedBase<T> &mat,
					  Eigen::Index outerIndex)
{
        std::vector<Eigen::Index> out;
	for (typename T::InnerIterator it(mat, outerIndex); it; ++it) {
	        out.push_back(T::IsRowMajor ? it.col() : it.row());
	}
	return out;
}

/* Since Eigen 3.4 is not supported in Rcpp, write own iterator */
template<typename T>
Eigen::MatrixXd matrix_slice_rows(const Eigen::DenseBase<T> &mat,
				  const std::vector<Eigen::Index> &indices)
{
        Eigen::MatrixXd out(indices.size(), mat.cols());
	for (Eigen::Index i : indices) {
	        out << mat.row(i).eval();
	}
	return out;
}
#endif
