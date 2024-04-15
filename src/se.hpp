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

/*
 * OLS Beta:
 * Use \beta = (X'X + \lambdaI)^{-1}X'Y
 *           = (V\Sigma'U'U\SigmaV' + \lambda I)^{-1}V\Sigma'U'U\SigmaV_i
 *           = (\V diag(\sigma^2+lambda)^{-1} V')V\Sigma^2V_i
 *           = \V diag(\sigma^2/(\sigma^2 + \lambda)) V_i
 *
 * When removing i from X,
 * (V_{-i} \Sigma U'U \Sigma V_{-i}' + \lambda I)^{-1} V_{-i}\Sigma'U'U\Sigma V_i
 * (V_{-i} (\Sigma^2 + \lambda) V_{-i}^{-1} V_{-i}\Sigma^2 V_i
 *
 * todo: use Sherman-Morrison: A^{-1} + A^{-1}vv'A^{-1}/(1 - v'A^{-1}v)
 *    = V\Sigma^{-2}V + (V\Sigma^{-2}V v v' V\Sigma^{-2}V)/(1 - v' V\Sigma^{-2}V v')
 *    = 
 */
template<typename TV>
Eigen::Vector<typename TV::Scalar, TV::ColsAtCompileTime>
intra_ols_beta(const Eigen::MatrixBase<TV> &V,
               const Eigen::DiagonalMatrix<typename TV::Scalar, TV::RowsAtCompileTime> &tuuituu,
               Eigen::Index i)
{
        Eigen::Vector<typename TV::Scalar, TV::ColsAtCompileTime> beta = V.transpose() * tuuituu * V.col(i);
        beta[i] = 0.; // Set to zero to make equivalent to subsetting V, but keeping dim.
        return beta;
}

/*
 * OLS Residual
 * Use residual = U * (V_i - V \beta)
 *
 * Convert to scalar of U for compatibility
 */
template<typename TU, typename TV>
Eigen::Vector<typename TU::Scalar, TU::RowsAtCompileTime>
intra_ols_resid(const Eigen::MatrixBase<TU> &U,
                const Eigen::MatrixBase<TV> &V,
                const Eigen::Vector<typename TV::Scalar, TV::ColsAtCompileTime> &beta,
                Eigen::Index i)
{
        Eigen::Vector<typename TU::Scalar, TV::RowsAtCompileTime> res_pc = (V.col(i) - V * beta).eval().template cast<typename TU::Scalar>();
        return U * res_pc;
}
/*
 */
template<typename TU>
Eigen::Matrix<float, TU::ColsAtCompileTime, TU::ColsAtCompileTime>
intra_ols_norm_mat(const Eigen::DiagonalMatrix<float, TU::ColsAtCompileTime> &tuui,
                   const Eigen::MatrixBase<TU> &U,
                   const Eigen::Vector<typename TU::Scalar, TU::RowsAtCompileTime> &resid,
                   Eigen::Matrix<float, TU::RowsAtCompileTime, TU::ColsAtCompileTime> &qr_buffer)
{
        // Using float for QR for caching/speed + space
        qr_buffer = (resid.asDiagonal() * U).template cast<float>().eval();
        // Use in-place QR to speed up
        Eigen::HouseholderQR<Eigen::Ref<Eigen::Matrix<float, TU::RowsAtCompileTime, TU::ColsAtCompileTime> > > qr(qr_buffer);
        Eigen::Matrix<float, TU::ColsAtCompileTime, TU::ColsAtCompileTime> R;
        R = qr.matrixQR().topLeftCorner(U.cols(), U.cols()).template triangularView<Eigen::Upper>();
        return R * tuui;
}

template<typename TV>
Eigen::Array<typename TV::Scalar, TV::ColsAtCompileTime, 1>
intra_ols_se(const Eigen::MatrixBase<TV> &V,
             const Eigen::Matrix<float, TV::RowsAtCompileTime, TV::RowsAtCompileTime> &shared,
             float epsilon=std::numeric_limits<float>::epsilon())
{
        // shared should be same dimension of R from qr(diag(resid) * U)
        Eigen::Matrix<typename TV::Scalar, TV::ColsAtCompileTime, TV::ColsAtCompileTime> shared_casted = shared.template cast<typename TV::Scalar>();
        Eigen::Vector<typename TV::Scalar, TV::ColsAtCompileTime> se = (shared_casted * V).colwise().norm().eval();
        return se.cwiseMax(epsilon).array();
}
template<typename TV>
std::pair<float, std::vector<Eigen::Index> > intra_ols_lambda(const Eigen::Array<float, TV::RowsAtCompileTime, 1> sigma,
                                              const Eigen::MatrixBase<TV> &V,
                                              Eigen::Index i,
                                              double edof,
                                              const Eigen::Array<float, TV::ColsAtCompileTime, 1> V_colnorm,
                                              float min_cor=0,
                                              bool abs_cutoff=false)
{
        Eigen::Vector<typename TV::Scalar, TV::ColsAtCompileTime> beta = V.transpose() * V.col(i);
        beta[i] = 0;
        Eigen::Vector<float, TV::RowsAtCompileTime> res_pc = (V.col(i) - V * beta).eval().template cast<float>();
        Eigen::Array<float, TV::RowsAtCompileTime, 1> res_pc_S = sigma * res_pc.array();
        double beta_norm2 = beta.squaredNorm();
        double lambda = res_pc_S.matrix().squaredNorm() / std::max(beta_norm2, 1e-100);
        lambda *= (V.cols() - 1.) / edof;
        std::vector<Eigen::Index> good;
        for (Eigen::Index j = 0; j < V.cols(); j++) {
                if (abs_cutoff && (beta(j) < V_colnorm(j) * V_colnorm(i) * min_cor * -1)) {
                        good.push_back(j);
                } else if (beta(j) > V_colnorm(j) * V_colnorm(i) * min_cor) {
                        good.push_back(j);
                }
        }
        return std::make_pair(static_cast<float>(lambda), good);
}

template<typename TU, typename TV>
Eigen::SparseMatrix<float> intra_robust_se(const Eigen::MatrixBase<TU> &U,
                                           const Eigen::MatrixBase<TV> &V,
                                           float min_cor=0,
                                           float t_cutoff=6.5,
                                           bool abs_cutoff=false)
{
        const Eigen::Array<float, TU::ColsAtCompileTime, 1> sigma_2 = U.colwise().squaredNorm().template cast<float>().eval().array();
        const Eigen::Array<float, TU::ColsAtCompileTime, 1> sigma_2_inv = sigma_2.inverse();
        const Eigen::DiagonalMatrix<float, TU::ColsAtCompileTime> tuui = sigma_2_inv.matrix().asDiagonal();

        Eigen::SparseMatrix<float> M(V.cols(), V.cols());
        const Eigen::Array<float, TV::ColsAtCompileTime, 1> V_colnorm = V.colwise().norm().array().template cast<float>().eval();
        ImplProgress p(V.cols());
#if defined(_OPENMP)
#pragma omp parallel 
#endif
        {
                Eigen::SparseMatrix<float> local_mat(V.cols(), V.cols());
                Eigen::Matrix<float, TU::RowsAtCompileTime, TU::ColsAtCompileTime> buffer;
#if defined(_OPENMP)
#pragma omp for nowait
#endif
                for (Eigen::Index i = 0; i < V.cols(); i++) {
                        if (!p.check_abort()) {
                                p.increment();
                                std::pair<float, std::vector<Eigen::Index> > pr = intra_ols_lambda(sigma_2, V, i, U.rows() - U.cols(), V_colnorm, min_cor, abs_cutoff);
                                float i_lambda = pr.first;
                                std::vector<Eigen::Index> good = pr.second;
                                if (good.size() == 0) { continue; }
                                const Eigen::DiagonalMatrix<typename TV::Scalar, TU::ColsAtCompileTime> tuuituu =
                                  (sigma_2 / (sigma_2 + i_lambda)).matrix().template cast<typename TV::Scalar>().asDiagonal();
                                Eigen::Vector<typename TV::Scalar, TV::ColsAtCompileTime> beta = intra_ols_beta(V, tuuituu, i);
                                Eigen::Matrix<float, TU::ColsAtCompileTime, TU::ColsAtCompileTime> shared = intra_ols_norm_mat(tuui, U, intra_ols_resid(U, V, beta, i), buffer);
                                Eigen::Array<typename TV::Scalar, TV::ColsAtCompileTime, 1> tv, se;
                                if (good.size()==V.cols()) {
                                        se = intra_ols_se(V, shared);
                                        tv = beta.array() / se;
                                } else {
                                        se = intra_ols_se(V(Eigen::all, good), shared);
                                        tv = beta(good).array() / se;
                                }
                                for (int k = 0; k < good.size(); k++) {
                                        Eigen::Index j = good[k];
                                        if (abs_cutoff && (t_cutoff <= -tv[k])) {
                                                local_mat.insert(j, i) = tv[k];
                                        } else if (tv[k] >= t_cutoff) {
                                                local_mat.insert(j, i) = tv[k];
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
