#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "se.hpp"

namespace py = pybind11;

template<typename TU, typename TV>
Eigen::SparseMatrix<float> py_robust_se(const py::EigenDRef<TU> &U,
					const py::EigenDRef<TV> &V,
					float lambda,
					float t_cutoff,
					bool abs_cutoff) 
{
	return intra_robust_se(U, V, lambda, t_cutoff, abs_cutoff);
}

PYBIND11_MODULE(_core, m) {
	m.def("py_robust_se_ff", &py_robust_se<Eigen::MatrixXf, Eigen::MatrixXf>);
	m.def("py_robust_se_df", &py_robust_se<Eigen::MatrixXd, Eigen::MatrixXf>);
	m.def("py_robust_se_fd", &py_robust_se<Eigen::MatrixXf, Eigen::MatrixXd>);
	m.def("py_robust_se_dd", &py_robust_se<Eigen::MatrixXd, Eigen::MatrixXd>);
}
