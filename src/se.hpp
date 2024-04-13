#ifndef SE_HPP
#define SE_HPP
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>


#include <algorithm>
#include <random>
#include <chrono>
#include <numeric>
#include <limits>
#include <vector>
#include "implprogress.hpp"

template<typename TV>
Eigen::Vector<typename TV::Scalar, Eigen::Dynamic>
intra_ols_beta(const Eigen::MatrixBase<TV> &V,
	       const Eigen::DiagonalMatrix<typename TV::Scalar, Eigen::Dynamic> &tuuituu,
	       Eigen::Index i)
{
	Eigen::Vector<typename TV::Scalar, Eigen::Dynamic> beta = (V.transpose() * tuuituu * V.col(i));
	beta[i] = 0.; // Set to zero to make equivalent to subsetting V, but keeping dim.
	return beta;
}

template<typename TU, typename TV>
Eigen::Vector<typename TU::Scalar, Eigen::Dynamic>
intra_ols_resid(const Eigen::MatrixBase<TU> &U,
		const Eigen::MatrixBase<TV> &V,
		const Eigen::Vector<typename TV::Scalar, Eigen::Dynamic> &beta,
		Eigen::Index i)
{
	Eigen::Vector<typename TU::Scalar, Eigen::Dynamic> res_pc = (V.col(i) - V * beta).eval().template cast<typename TU::Scalar>();
	return U * res_pc;
}
template<typename TU>
Eigen::MatrixXf intra_ols_norm_mat(const Eigen::DiagonalMatrix<float, Eigen::Dynamic> &tuui,
				   const Eigen::MatrixBase<TU> &U,
				   const Eigen::Vector<typename TU::Scalar, Eigen::Dynamic> &resid)
{
	// Using float for QR for caching/speed + space
	Eigen::Index n = U.cols();
	Eigen::MatrixXf pre_qr = (resid.asDiagonal() * U).template cast<float>().eval();
	// Use in-place QR to speed up
	Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXf> > qr(pre_qr);
	Eigen::MatrixXf R = qr.matrixQR().topLeftCorner(n, n).triangularView<Eigen::Upper>();
	return R * tuui;
}

template<typename TV>
Eigen::Array<typename TV::Scalar, Eigen::Dynamic, 1> intra_ols_se(const Eigen::MatrixBase<TV> &V,
								  const Eigen::MatrixXf &shared,
								  float epsilon=std::numeric_limits<float>::epsilon())
{
	// shared should be same dimension of R from qr(diag(resid) * U)
	Eigen::Matrix<typename TV::Scalar, Eigen::Dynamic, Eigen::Dynamic> shared_casted = shared.template cast<typename TV::Scalar>();
	Eigen::Vector<typename TV::Scalar, Eigen::Dynamic> se = (shared_casted * V).colwise().norm().eval();
	return se.cwiseMax(epsilon).array();
}
template<typename TU, typename TV>
Eigen::SparseMatrix<float> intra_robust_se(const Eigen::MatrixBase<TU> &U,
					   const Eigen::MatrixBase<TV> &V,
					   float lambda=100,
					   float t_cutoff=6.5,
					   bool abs_cutoff=false)
{
	const Eigen::VectorXf sigma_2 = U.colwise().squaredNorm().template cast<float>().eval();
	const Eigen::ArrayXf sigma_2_inv = (1. / (sigma_2.array() + lambda*lambda));
	const Eigen::DiagonalMatrix<float, Eigen::Dynamic> tuui = sigma_2_inv.matrix().asDiagonal();
	const Eigen::DiagonalMatrix<typename TV::Scalar, Eigen::Dynamic> tuuituu = (sigma_2.array() / (sigma_2.array() + lambda*lambda)).matrix().template cast<typename TV::Scalar>().asDiagonal();
	Eigen::SparseMatrix<float> M(V.cols(), V.cols());
	ImplProgress p(V.cols());
#if defined(_OPENMP)
#pragma omp parallel 
#endif
	{
		Eigen::SparseMatrix<float> local_mat(V.cols(), V.cols());
#if defined(_OPENMP)
#pragma omp for nowait
#endif
		for (Eigen::Index i = 0; i < V.cols(); i++) {
			if (!p.check_abort()) {
				p.increment();
				Eigen::Vector<typename TV::Scalar, Eigen::Dynamic> beta = intra_ols_beta(V, tuuituu, i);
				Eigen::MatrixXf shared = intra_ols_norm_mat(tuui, U, intra_ols_resid(U, V, beta, i));
				Eigen::Array<typename TV::Scalar, Eigen::Dynamic, 1> se = intra_ols_se(V, shared);
				Eigen::Array<typename TV::Scalar, Eigen::Dynamic, 1> tv = beta.array() / se;
				for (int j = 0; j < tv.size(); j++) {
					if (abs_cutoff && (t_cutoff <= -tv[j])) {
						local_mat.insert(j, i) = tv[j];
					} else if (tv[j] >= t_cutoff) {
						local_mat.insert(j, i) = tv[j];
					}
				}
			}
		}
#if defined(_OPENMP)
#pragma omp critical
#endif
		M += local_mat;
	}
	if (!p.check_abort()) {
		M.makeCompressed();
		return M;
	} else {
		Eigen::SparseMatrix<float> Mbad(V.cols(), V.cols());
		return Mbad;
	}
}

#endif
