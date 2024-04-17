#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "se.hpp"

namespace py = pybind11;

template<typename TU, typename TV>
Eigen::SparseMatrix<float> py_robust_se(const py::EigenDRef<TU> &U,
					const py::EigenDRef<TV> &V,
					float min_cor,
					float t_cutoff,
					bool abs_cutoff,
					float lambda_pow) 
{
	return intra_robust_se(U, V, min_cor, t_cutoff, abs_cutoff, lambda_pow);
}

PYBIND11_MODULE(_core, m) {
	m.def("py_robust_se_ff", &py_robust_se<Eigen::MatrixXf, Eigen::MatrixXf>);
	m.def("py_robust_se_df", &py_robust_se<Eigen::MatrixXd, Eigen::MatrixXf>);
	m.def("py_robust_se_fd", &py_robust_se<Eigen::MatrixXf, Eigen::MatrixXd>);
	m.def("py_robust_se_dd", &py_robust_se<Eigen::MatrixXd, Eigen::MatrixXd>);
}
