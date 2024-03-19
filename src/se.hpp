#ifndef SE_HPP
#define SE_HPP
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

//#if defined(_USE_GSL)
// #include <gsl/gsl_statistics.h>
// #include <gsl/gsl_cdf.h>
//#endif

#include <algorithm>
#include <random>
#include <chrono>
#include <numeric>
#include <vector>
#include "implprogress.hpp"
/*
 * Reminder: Y has many columns, X should have 1.
 */
template<typename TX, typename TY>
Eigen::MatrixXd ols_beta(const Eigen::MatrixBase<TX> &X, const Eigen::MatrixBase<TY> &Y, double lambda)
{
  	if (lambda > 0) {
        	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(X.cols(), X.cols());
                Eigen::MatrixXd A = (X.transpose() * X + lambda * I).eval();
                Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
                return svd.solve(X.transpose() * Y);
        } else {
        	return X.completeOrthogonalDecomposition().solve(Y).eval();
        }
}

template<typename TX, typename TY, typename TB>
Eigen::MatrixXd ols_resid(const Eigen::MatrixBase<TX> &X,
			  const Eigen::MatrixBase<TY> &Y,
			  const Eigen::MatrixBase<TB> &beta) {
	return (Y - X * beta).eval();
}

/*
 * Y must be U1^+ U1) 
 */
template<typename TX, typename TY>
Eigen::ArrayXd robust_se_X(const Eigen::MatrixBase<TX> &Xmat,
			   const Eigen::MatrixBase<TY> &Y,
			   double lambda=1e-10,
			   double epsilon=1e-300)
{
	Eigen::VectorXd X = Xmat.col(0).eval();
	double X_sqnorm = std::max(X.squaredNorm(), epsilon);
	// X pseudoinverse is X' / squared_norm
	Eigen::VectorXd beta = (X.transpose() * Y).eval() / (X_sqnorm + lambda);
	Eigen::VectorXd var = ((X.transpose() / X_sqnorm).cwiseAbs2() * (Y - X * beta.transpose()).cwiseAbs2()).eval();
	Eigen::ArrayXd tval = beta.array();
	return tval * var.cwiseMax(epsilon).array().rsqrt();
}

#include <iostream>
/* Each column of X versus Y */
template<typename TX, typename TY>
Eigen::ArrayXd robust_se_Y(const Eigen::MatrixBase<TX> &X,
			   const Eigen::MatrixBase<TY> &Ymat,
			   double lambda=1e-10,
			   double epsilon=1e-300)
{
  
	Eigen::VectorXd Y = Ymat.col(0);
        Eigen::ArrayXd X_sqnorm = X.colwise().squaredNorm().cwiseMax(epsilon).array();
	// X pseudoinverse is X' / squared_norm
	Eigen::VectorXd beta = ((X.transpose() * Y).array() / (X_sqnorm + lambda)).matrix().eval();
	Eigen::MatrixXd res2 = X * beta.asDiagonal();
	res2.colwise() -= Y;
	res2 = res2.cwiseAbs2();
	//Eigen::MatrixXd res2 = (Y - (X * beta.asDiagonal()).colwise()).cwiseAbs2().eval();
	std::cout << "res2 R:" << res2.rows() << " C:" << res2.cols()  << std::endl;
	Eigen::MatrixXd pinv2 = (X * (1/X_sqnorm).matrix().asDiagonal()).cwiseAbs2().eval();
	std::cout << "pinv2 R:" << pinv2.rows() << " C:" << pinv2.cols()  << std::endl;
	Eigen::VectorXd var = (pinv2.array() * res2.array()).colwise().sum();
	std::cout << "var size:" << var.size() << std::endl;
	//Eigen::VectorXd var = ((X * (1/X_sqnorm).matrix().asDiagonal()).transpose().cwiseAbs2() * (Y - X * beta.asDiagonal()).cwiseAbs2()).eval();
	Eigen::ArrayXd tval = beta.array();
	return tval * var.cwiseMax(epsilon).array().rsqrt();
}
template<typename TY, typename TL>
Eigen::VectorXd ols_beta_L(Eigen::VectorXd X, double X_sqnorm,
                           const Eigen::MatrixBase<TY> &Y,
                           const Eigen::ArrayBase<TL> &lambda)
{
  	Eigen::ArrayXd numer = (X.transpose() * Y).eval().array();
        Eigen::ArrayXd denom = (lambda + X_sqnorm).eval();
        Eigen::ArrayXd quot = numer / denom;
        return quot.matrix();
}
template<typename TX, typename TY, typename TL>
Eigen::ArrayXd robust_se_L(const Eigen::MatrixBase<TX> &Xmat,
			   const Eigen::MatrixBase<TY> &Y, /* V\Sigma?? */
                           const Eigen::ArrayBase<TL> &lambda,
			   double epsilon=1e-300)
{
	Eigen::VectorXd X = Xmat.col(0).eval();
	double X_sqnorm = std::max(X.squaredNorm(), epsilon);
        Eigen::VectorXd beta = ols_beta_L(X, X_sqnorm, Y, lambda);
	// X pseudoinverse is X' / squared_norm
	Eigen::VectorXd var = ((X.transpose() / X_sqnorm).cwiseAbs2() * (Y - X * beta.transpose()).cwiseAbs2()).eval();
	Eigen::ArrayXd tval = beta.array();
	return tval * var.cwiseMax(epsilon).array().rsqrt();
}
template<typename TX, typename TY>
Eigen::SparseMatrix<double> robust_se(const Eigen::MatrixBase<TX> &X,
				      const Eigen::MatrixBase<TY> &Y,
				      double lambda=1e-10,
				      double epsilon=1e-300,
				      double t_cutoff=6.5,
				      bool abs_cutoff=false)
{
	Eigen::SparseMatrix<double> M(Y.cols(), X.cols());
	ImplProgress p(X.cols());
#if defined(_OPENMP)
#pragma omp parallel 
#endif
	{
		Eigen::SparseMatrix<double> local_mat(Y.cols(), X.cols());
#if defined(_OPENMP)
#pragma omp for nowait
#endif
		for (int i = 0; i < X.cols(); i++) {
			if (!p.check_abort()) {
		        	p.increment();
				Eigen::ArrayXd tv = robust_se_X(X.col(i), Y, lambda, epsilon);
				for (int j = 0; j < tv.size(); j++) {
					if (abs_cutoff && (t_cutoff <= -tv[j])) {
						local_mat.insert(j, i) = -tv[j];
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
		Eigen::SparseMatrix<double> Mbad(Y.cols(), X.cols());
		return Mbad;
	}
}

#endif
