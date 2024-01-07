#ifndef SE_HPP
#define SE_HPP
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <algorithm>

/*
 * Reminder: Y has many columns, X should have 1.
 */
template<typename TX, typename TY>
Eigen::MatrixXd ols_beta(const Eigen::MatrixBase<TX> &X, const Eigen::MatrixBase<TY> &Y)
{
    return X.completeOrthogonalDecomposition().solve(Y).eval();
}

template<typename TX, typename TY, typename TB>
Eigen::MatrixXd ols_resid(const Eigen::MatrixBase<TX> &X,
			  const Eigen::MatrixBase<TY> &Y,
			  const Eigen::MatrixBase<TB> &beta) {
	return (Y - X * beta).eval();
}

/*
 * Frisch-Waugh partialling out.
 * Since I - B(B'B)B' is multiplied by U,
 * U - B'(B'B)B'U will be the new transform.
 * = U - B (B^{+} U)
 * = U - B BPU (where BPU=B^{+} U, saving nrow(U) dimensions ).
 * Additionally, allows BPU to be calculated using COD/pivoted Householder QR which will be more stable than B^{+} * U
 */
template<typename TU, typename TB, typename TBPU>
Eigen::VectorXd adjust_res(const Eigen::MatrixXd &res,
			   const Eigen::MatrixBase<TU> &U,
			   const Eigen::MatrixBase<TB> &B,
			   const Eigen::MatrixBase<TBPU> &BPU) {
	/* First coerce bpu_res as ncol(B) is likely very small */
	Eigen::MatrixXd bpu_res = (BPU * res).eval();
	return (U * res - B * bpu_res).colwise().norm().eval();
}
/*
 * Ideally this would be called on a KNN graph, such that each X has K neighbors in Vs.
 * BPU should be computed in the outer layer
 */
template<typename TX, typename TVs, typename TU, typename TTUU, typename TB, typename TBPU>
Eigen::VectorXd hc0_se_Xvec(const Eigen::MatrixBase<TX> &X,
			    const Eigen::MatrixBase<TVs> &Vs, /* Vs */
			    const Eigen::MatrixBase<TU> &U,
			    const Eigen::MatrixBase<TTUU> &TUU,
			    const Eigen::MatrixBase<TB> &B,
			    const Eigen::MatrixBase<TBPU> &BPU) {
	const Eigen::MatrixXd betas = ols_beta(X, Vs);
	Eigen::MatrixXd res = ols_resid(X, Vs, betas);
	Eigen::MatrixXd bread_inv = (X.transpose() * TUU * X).eval();
	Eigen::VectorXd ssr = adjust_res(res, U, B, BPU);
	Eigen::VectorXd se = ssr.cwiseSqrt() / bread_inv(0, 0);
	return betas.row(0).transpose().cwiseQuotient(se).eval();
}
#endif
