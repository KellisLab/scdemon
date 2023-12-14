
#include <Eigen/Core>
#include <Eigen/Dense>

Eigen::MatrixXd ols(const Eigen::VectorXd &X, const Eigen::MatrixXd &Y)
{
	return (X.transpose() * X).inverse() * X.transpose() * Y;
}
