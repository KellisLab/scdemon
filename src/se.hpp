#ifndef SE_HPP
#define SE_HPP
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

//#if defined(_USE_GSL)
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cdf.h>
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

template<typename T>
std::vector<T> sample(int N)
{
	std::vector<T> seq(N);
        std::iota(seq.begin(), seq.end(), 1);
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle(seq.begin(), seq.end(), std::default_random_engine(seed));
        return seq;
}
#include <iostream>
template<typename TU, typename TV>
Eigen::ArrayXd robust_se_Xfull(const Eigen::Index &x_idx,
			       const Eigen::MatrixBase<TU> &U,
			       const Eigen::MatrixBase<TV> &V,
			       double epsilon=1e-300,
			       Eigen::Index block_size=10000)
{
	Eigen::VectorXd UX = (U * V.col(x_idx)).eval();
        double UX_rsqnorm = 1/std::max(UX.squaredNorm(), epsilon);
	Eigen::RowVectorXd UXp = (UX.transpose() * UX_rsqnorm);
        Eigen::RowVectorXd UXp2 = UXp.cwiseAbs2();
        Eigen::RowVectorXd beta = (UXp * U) * V;
        Eigen::MatrixXd meat_nou = V - V.col(x_idx) * beta;
        Eigen::RowVectorXd var = Eigen::RowVectorXd::Zero(V.cols());
        std::vector<Eigen::Index> U_idx = sample<Eigen::Index>(U.rows());
        Eigen::Index n_sample = ceil(sqrt(U.rows()));
        Eigen::Index nblocks = n_sample / block_size + (n_sample % block_size != 0);
        for (Eigen::Index i = 0; i < nblocks; i++) {
		Eigen::Index left = i * block_size;
		Eigen::Index right = std::min(n_sample, (i+1)*block_size);
                Eigen::MatrixXd Ui(right - left, U.cols());
                Eigen::RowVectorXd UXp2i = Eigen::RowVectorXd::Zero(right - left);
                for (Eigen::Index j = 0; j < Ui.rows(); j++) {
                	UXp2i(j) = UXp2(U_idx[j + left]);
                        Ui.row(j) = U.row(U_idx[j + left]);
                }
        	var += UXp2i * (Ui * meat_nou).cwiseAbs2();
        }
        var *= U.rows() / n_sample; // scaling because sampled fraction of U.rows()
        return beta.array() * var.cwiseMax(epsilon).array().rsqrt();
}
template<typename TU, typename TV>
Eigen::ArrayXd robust_se_Xfull_bad(const Eigen::Index &x_idx,
			       const Eigen::MatrixBase<TU> &U,
			       const Eigen::MatrixBase<TV> &V,
			       double epsilon=1e-300,
			       Eigen::Index block_size=100)
{
	Eigen::VectorXd X = (U * V.col(x_idx)).eval();
	double X_rsqnorm = 1/std::max(X.squaredNorm(), epsilon);
	Eigen::VectorXd Xp2T = (X * X_rsqnorm).cwiseAbs2();
        // take nblocks as ceiling of ncol / block_size
        Eigen::Index nblocks = V.cols() / block_size + (V.cols() % block_size != 0);
        Eigen::ArrayXd var(V.cols());
	for (Eigen::Index i = 0; i < nblocks; i++) {
        	std::cout << i << std::endl;
		Eigen::Index left = i * block_size;
		Eigen::Index right = std::min(V.cols(), (i+1)*block_size);
                std::vector<Eigen::Index> cols(right - left);
                std::iota(cols.begin(), cols.end(), left);
                Eigen::MatrixXd V1(V.rows(), cols.size());
                for (Eigen::Index j = 0; j < cols.size(); j++) {
                	V1.col(j) = V.col(cols[j]);
                }
		Eigen::MatrixXd Y = U * V1;
		Y = (Y - X * X.transpose() * Y * X_rsqnorm).cwiseAbs2();
                Eigen::VectorXd bvar = Xp2T.transpose() * Y;
                for (Eigen::Index j = 0; j < cols.size(); j++) {
                	var(cols[j]) = bvar(j);
                }
	}
        Eigen::ArrayXd tval = X.array() * X_rsqnorm * var.cwiseMax(epsilon).array().rsqrt();
        return tval;
}
/*
 * Y must be U1^+ U1) 
 */
template<typename TY>
Eigen::ArrayXd robust_se_X(const Eigen::Index &x_idx,
			   const Eigen::MatrixBase<TY> &Y, /* V\Sigma?? */
			   double epsilon=1e-300)
{
  // dimensions
	Eigen::VectorXd X = Y.col(x_idx).eval();
	double X_rsqnorm = 1/std::max(X.squaredNorm(), epsilon);
	// X pseudoinverse is X' / squared_norm
	Eigen::VectorXd beta = (X.transpose() * Y).eval() * X_rsqnorm;
	Eigen::VectorXd var = ((X.transpose() * X_rsqnorm).cwiseAbs2() * (Y - X * beta.transpose()).cwiseAbs2()).eval();
	Eigen::ArrayXd tval = beta.array();
	return tval * var.cwiseMax(epsilon).array().rsqrt();
}

template<typename TY>
Eigen::SparseMatrix<double> robust_se(const Eigen::MatrixBase<TY> &Y,
				      implprogress_callback callback,
				      implprogress_callback interrupt_checker,
				      double epsilon=1e-300,
				      double t_cutoff=6.5,
				      bool abs_cutoff=false)
{
	Eigen::SparseMatrix<double> M(Y.cols(), Y.cols());
	ImplProgress p(Y.cols(), callback, interrupt_checker);
#if defined(_OPENMP)
#pragma omp parallel 
#endif
	{
		Eigen::SparseMatrix<double> local_mat(Y.cols(), Y.cols());
#if defined(_OPENMP)
#pragma omp for nowait
#endif
		for (int i = 0; i < Y.cols(); i++) {
			if (!p.check_abort()) {
		        	p.increment();
				Eigen::ArrayXd tv = robust_se_X(i, Y, epsilon);
				for (int j = 0; j < tv.size(); j++) {
					if (i != j) {
						if (abs_cutoff && (t_cutoff <= -tv[j])) {
							local_mat.insert(j, i) = -tv[j];
						} else if (tv[j] >= t_cutoff) {
							local_mat.insert(j, i) = tv[j];
						}
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
		Eigen::SparseMatrix<double> Mbad(Y.cols(), Y.cols());
		return Mbad;
	}
}

template<typename T>
Eigen::ArrayXd cwiseVar(const Eigen::MatrixBase<T> &Y)
{
	Eigen::VectorXd means = Y.colwise().mean().eval();
	return (Y - means).array().square().colwise().sum() / Y.rows();
}

template<typename TY, typename TD>
Eigen::SparseMatrix<double> robust_se_pvalue(const Eigen::MatrixBase<TY> &Y,
					     const Eigen::ArrayBase<TD> &dof,
					     implprogress_callback callback,
					     implprogress_callback interrupt_checker,
					     double nominal_p_cutoff=0.05,
					     bool abs_cutoff=false,
					     double epsilon=1e-300)
{
	Eigen::SparseMatrix<double> M(Y.cols(), Y.cols());
	ImplProgress p(Y.cols(), callback, interrupt_checker);
	double adj_p_cutoff = nominal_p_cutoff / (Y.cols()*Y.cols());
	Eigen::ArrayXd t_cutoff(Y.cols());
	for (int i = 0; i < Y.cols(); i++) {
		t_cutoff(i) = gsl_cdf_tdist_Qinv(adj_p_cutoff, dof(i));
	}
#if defined(_OPENMP)
#pragma omp parallel 
#endif
	{
		Eigen::SparseMatrix<double> local_mat(Y.cols(), Y.cols());
#if defined(_OPENMP)
#pragma omp for nowait
#endif
		for (int i = 0; i < Y.cols(); i++) {
			if (!p.check_abort()) {
				p.increment();
				Eigen::ArrayXd tv = robust_se_X(i, Y, epsilon);
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
							local_mat.insert(j, i) = std::max(pval, epsilon);
						}
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
		Eigen::SparseMatrix<double> Mbad(Y.cols(), Y.cols());
		return Mbad;
	}
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
