#ifndef SE_HPP
#define SE_HPP
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cdf.h>

#include <algorithm>
#include <vector>

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
Eigen::VectorXd fw_meat(const Eigen::MatrixXd &res,
			const Eigen::MatrixBase<TU> &U,
			const Eigen::MatrixBase<TB> &B,
			const Eigen::MatrixBase<TBPU> &BPU) {
	/* First coerce bpu_res as ncol(B) is likely very small */
  /* todo figure out FW bread that uses the transform I - X(XX)X here */
  /* Can we pre-multiply res times X? */
	Eigen::MatrixXd bpu_res = (BPU * res).eval();
	return (U * res - B * bpu_res).colwise().squaredNorm().eval();
}

template<typename TU, typename TB, typename TBPU>
long double fw_bread(const Eigen::MatrixXd &X,
		     const Eigen::MatrixBase<TU> &U,
		     const Eigen::MatrixBase<TB> &B,
		     const Eigen::MatrixBase<TBPU> &BPU) {
	/* First coerce bpu_res as ncol(B) is likely very small */
  /* todo figure out FW bread that uses the transform I - X(XX)X here */
  /* bread: figure out once: U * X - B * bpu_X*/
	Eigen::MatrixXd bpu_X = (BPU * X).eval();
	Eigen::MatrixXd adjX = U * X - B * bpu_X;
	Eigen::MatrixXd xtx = adjX.transpose() * adjX;
	return 1 / xtx(0, 0);
}
/*
 * Ideally this would be called on a KNN graph, such that each X has K neighbors in Vs.
 * BPU should be computed in the outer layer.
 * TODO use covariate free model for KNN graph?
 * TODO maybe use samples of U;B (but not BPU) to calculate initial graph
 * TODO blocks of X as well?
 * TODO maybe simplify computation: (X^+ Y) / 
 */
template<typename TX, typename TVs, typename TU, typename TB, typename TBPU>
Eigen::VectorXd old_hc0_se_Xvec(const Eigen::MatrixBase<TX> &X,
			    const Eigen::MatrixBase<TVs> &Vs, /* Vs */
			    const Eigen::MatrixBase<TU> &U,
			    const Eigen::MatrixBase<TB> &B,
			    const Eigen::MatrixBase<TBPU> &BPU,
			    Eigen::Index blockSize = 2048) {
  /* Todo: for future,
   * have just Vs, X_ind, and Y_ind,
   * and calculate betas, resid iteratively.
   * Here, take row 0 because X should be 1 feature.
   */
	const Eigen::VectorXd betas = ols_beta(X, Vs).row(0); 
	/* Since X should be a column vector, res should be a matrix as well,
	 * of dimension {length(X), length(Vs)}, i.e. 50 x 36k (ncol Vs)
	 * TODO correct resid with dof n/(n-k)
	 */
	Eigen::MatrixXd resid = ols_resid(X, Vs, betas); /* parallelization: when more X, compute in block loop so still matrix not tensor */
	/* This assumes X has only one column, i.e. {U.nrow, 1} */
	Eigen::Matrix bpu_X = (BPU * X).eval();
	Eigen::VectorXd Xfw = (U * X - B * bpu_X).eval(); /* target for parallelization: more X */

	/* X_pinv_denom is the denominator of (Xfw+)**2,
	 * that is, Xfw+=Xfw.T / (Xfw.T * Xfw), and Xfw.T * Xfw is a scalar,
	 * so the squared scalar denom is computed here.
	 */
	double X_pinv_denom = Xfw.squaredNorm();
	bpu_X.resize(0, 0); /* free */
	/* Product of squares: Xfw^2 * (resid.fw)^2
	 * which is a reordering of Xfw+ * diag(resid.fw^2) * (Xfw+).T such that
	 * products are vectorized over all resid instead of per-resid
	 */
	Eigen::MatrixXd bpu_resid = (BPU * resid).eval();
	// Eigen::VectorXd var_times_denom = (Xfw.transpose().cwiseAbs2() * (U * resid - B * bpu_resid).cwiseAbs2()).eval();

	Eigen::VectorXd var_times_denom = Eigen::VectorXd::Zero(resid.cols());
	for (Eigen::Index left = 0; left < U.rows(); left += blockSize) {
	  Eigen::Index right = std::min(blockSize, U.rows() - left);
	  var_times_denom += Xfw.segment(left, right).transpose().cwiseAbs2() * (U.block(left, 0, right, U.cols()) * resid 
										 - B.block(left, 0, right, B.cols()) * bpu_resid).cwiseAbs2();
	}
	Eigen::VectorXd partial_se = var_times_denom.cwiseSqrt().cwiseMax(1e-300).eval();
	bpu_resid.resize(0, 0); /* free */
	/* Note that */
	return (betas * X_pinv_denom).cwiseQuotient(partial_se).eval();
}

template<typename TA, typename TB>
Eigen::MatrixXd cbind(const Eigen::MatrixBase<TA> &A,
		      const Eigen::MatrixBase<TB> &B)
{
	Eigen::MatrixXd C(A.rows(), A.cols() + B.cols());
	C.leftCols(A.cols()) = A;
	C.rightCols(B.cols()) = B;
	return C;
}
/*
 * TODO run regression on B first, then FW on residuals
 * OLS y ~ B 
 * so use B (BPU Vs) to get resid
 */
template<typename TY, typename TU, typename TB>
Eigen::ArrayXd robust_se_X(const Eigen::Index &x_idx,
			   const Eigen::MatrixBase<TY> &Y, /* V\Sigma?? */
			   const Eigen::MatrixBase<TU> &UpU,
			   const Eigen::MatrixBase<TB> &UpB,
			   double epsilon=1e-300)
{
	Eigen::MatrixXd UpBX = cbind(UpU * Y.col(x_idx), UpB);
	Eigen::MatrixXd BXpU = UpBX.completeOrthogonalDecomposition().pseudoInverse();
	Eigen::MatrixXd Ur_pre = UpU - UpBX * BXpU; // Annihilator matrix in U space
	Eigen::VectorXd var = BXpU.row(0).cwiseAbs2() * (Ur_pre * Y).cwiseAbs2();
	Eigen::ArrayXd tval = (BXpU.row(0) * Y).array();
	return tval * var.cwiseMax(epsilon).array().rsqrt();
}

template<typename TY, typename TU, typename TB>
Eigen::SparseMatrix<double> robust_se(const Eigen::MatrixBase<TY> &Y,
				      const Eigen::MatrixBase<TU> &UpU,
				      const Eigen::MatrixBase<TB> &UpB,
				      double epsilon=1e-300,
				      double t_cutoff=6.5,
				      bool abs_cutoff=false)
{
	Eigen::SparseMatrix<double> M(Y.cols(), Y.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 100)
#endif
	for (int i = 0; i < Y.cols(); i++) {
		Eigen::ArrayXd tv = robust_se_X(i, Y, UpU, UpB, epsilon);
		Eigen::SparseMatrix<double> MR(Y.cols(), Y.cols());
		for (int j = 0; j < tv.size(); j++) {
			if (i != j) {
				if (abs_cutoff && (t_cutoff <= -tv[j])) {
					MR.insert(j, i) = -tv[j];
				} else if (tv[j] >= t_cutoff) {
					MR.insert(j, i) = tv[j];
				}
			}
		}
#if defined(_OPENMP)
#pragma omp critical
#endif
		{
		  if (i % 100 == 0) { std::cout << i << "/" << Y.cols() << std::endl; }
		  M += MR;
		}
	}
	M.makeCompressed();
	return M;
}

template<typename T>
Eigen::ArrayXd cwiseVar(const Eigen::MatrixBase<T> &Y)
{
	Eigen::VectorXd means = Y.colwise().mean().eval();
	return (Y - means).array().square().colwise().sum() / Y.rows();
}

template<typename TY, typename TU, typename TB, typename TD>
Eigen::SparseMatrix<double> robust_se_pvalue(const Eigen::MatrixBase<TY> &Y,
					     const Eigen::MatrixBase<TU> &UpU,
					     const Eigen::MatrixBase<TB> &UpB,
					     const Eigen::ArrayBase<TD> &dof,
					     double nominal_p_cutoff=0.05,
					     bool abs_cutoff=false,
					     double epsilon=1e-300)
{
	Eigen::SparseMatrix<double> M(Y.cols(), Y.cols());
	// Eigen::ArrayXd var = cwiseVar(Y) / (dof + 1);
	double adj_p_cutoff = nominal_p_cutoff / (Y.cols()*Y.cols());
	Eigen::ArrayXd t_cutoff(Y.cols());
	for (int i = 0; i < Y.cols(); i++) {
		t_cutoff(i) = gsl_cdf_tdist_Qinv(adj_p_cutoff, dof(i));
	}
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 100)
#endif
	for (int i = 0; i < Y.cols(); i++) {
		// Welch–Satterthwaite equation
		//Eigen::ArrayXd adj_dof = (var + var(i)).square() / (var.square() / dof + var(i)*var(i)/dof(i)).max(epsilon);
		Eigen::ArrayXd tv = robust_se_X(i, Y, UpU, UpB, epsilon);
		Eigen::SparseMatrix<double> MR(Y.cols(), Y.cols());
		for (int j = 0; j < tv.size(); j++) {
			if (i != j) {
				double pval = 1;
				if (abs_cutoff && (tv[j] <= -t_cutoff(j))) {
					pval = gsl_cdf_tdist_P(tv[j], dof(j));
				} else if (tv[j] >= t_cutoff(j)) {
					pval = gsl_cdf_tdist_Q(tv[j], dof(j));
				}
				if (abs_cutoff) {
					pval = 2 * pval;
				}
				if (pval < adj_p_cutoff) {
					MR.insert(j, i) = std::max(pval, epsilon);
				}
			}
		}
#if defined(_OPENMP)
#pragma omp critical
#endif
		{
		  if (i % 100 == 0) { std::cout << i << "/" << Y.cols() << std::endl; }
		  M += MR;
		}
	}
	M.makeCompressed();
	return M;
}
/* Usage: Make sure to check if empty before slicing */
template<class T>
std::vector<Eigen::Index> sparse_extract_inner(const Eigen::SparseCompressedBase<T> &mat,
					       Eigen::Index outerIndex)
{
	std::vector<Eigen::Index> inner_nnz;
	if (T::IsRowMajor) {
		for (typename T::InnerIterator it(mat, outerIndex); it; ++it) {
			inner_nnz.push_back(it.col());
		}
	} else {
		for (typename T::InnerIterator it(mat, outerIndex); it; ++it) {
			inner_nnz.push_back(it.row());
		}
	}
	return inner_nnz;
}
#endif
