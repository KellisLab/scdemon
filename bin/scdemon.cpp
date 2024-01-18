
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "se.hpp"
int main(int argc, char**argv)
{
  //const char* const short_opts = "i:u:v:o:";
  // Eigen::MatrixXd X_pca = h5read(argv[1], "/obsm/X_pca");
  Eigen::MatrixXd X_pca = Eigen::MatrixXd::Random(50000, 50);
  Eigen::MatrixXd PCs = Eigen::MatrixXd::Random(50, 10000);
  Eigen::MatrixXd UpU = ols_beta(X_pca, X_pca);
  Eigen::MatrixXd B = Eigen::MatrixXd::Constant(50000, 1, 1.0);
  Eigen::MatrixXd UpB = ols_beta(X_pca, B);
  Eigen::ArrayXd dof = 50000 * Eigen::ArrayXd::Ones(10000);
  Eigen::SparseMatrix<double> P = robust_se_pvalue(PCs, UpU, UpB, dof, 0.05, false, 1e-300);
  // Eigen::VectorXd sigma = X_pca.colwise().norm();
  // Eigen::MatrixXd PCs = h5read(argv[1], "/varm/PCs").transpose();
  std::cout << P.rows() << " " << P.cols() << std::endl;
}
