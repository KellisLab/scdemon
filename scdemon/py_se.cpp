#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "se.hpp"

namespace py = pybind11;


// todo fix by using .cast<float>() in se.hpp
//template<typename TY, typename TU, typename TB>
Eigen::ArrayXd py_robust_se_X(int x_idx,
			      const py::EigenDRef<Eigen::MatrixXd> &Y,
			      const py::EigenDRef<Eigen::MatrixXd> &UpU,
			      const py::EigenDRef<Eigen::MatrixXd> &UpB)
{
	return robust_se_X(x_idx, Y, UpU, UpB, 1e-300);
}

Eigen::SparseMatrix<double> py_robust_se(const py::EigenDRef<Eigen::MatrixXd> &Y,
					 const py::EigenDRef<Eigen::MatrixXd> &UpU,
					 const py::EigenDRef<Eigen::MatrixXd> &UpB,
					 double t_cutoff=6.5,
					 bool abs_t=false)
{
	return robust_se(Y, UpU, UpB, 1e-300, t_cutoff, abs_t);
}
PYBIND11_MODULE(scdemon_ext, m) {
	//m.def("robust_se", &robust_se, "Robust SE");
	m.def("robust_se_X", &py_robust_se_X);
	m.def("robust_se", &py_robust_se);
	// <Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>);
	// m.def("robust_se_X_ddf", &py_robust_se_X<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXf>);
	// m.def("robust_se_X_dfd", &py_robust_se_X<Eigen::MatrixXd, Eigen::MatrixXf, Eigen::MatrixXd>);
	// m.def("robust_se_X_dff", &py_robust_se_X<Eigen::MatrixXd, Eigen::MatrixXf, Eigen::MatrixXf>);
	// m.def("robust_se_X_fdd", &py_robust_se_X<Eigen::MatrixXf, Eigen::MatrixXd, Eigen::MatrixXd>);
	// m.def("robust_se_X_fdf", &py_robust_se_X<Eigen::MatrixXf, Eigen::MatrixXd, Eigen::MatrixXf>);
	// m.def("robust_se_X_ffd", &py_robust_se_X<Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXd>);
	// m.def("robust_se_X_fff", &py_robust_se_X<Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf>);  
}
