
#include <Eigen/Core>
#include <Eigen/Dense>
#include <memory>
#include <algorithm>

template<typename T>
Eigen::MatrixXd ols_beta(const Eigen::VectorXd &X, const Eigen::MatrixXd &Y)
{
	return (X.transpose() * X).inverse() * X.transpose() * Y;
}

Eigen::VectorXd hc0_se_Xvec(const Eigen::VectorXd &X,
			    const Eigen::MatrixXd &Y,
			    const std::shared_ptr<Eigen::VectorXd> &dofPtr = nullptr,
			    const std::shared_ptr<Eigen::MatrixXd> &transform = nullptr)
{
	/*
	 * todo: change X_pinv to be in U space
	 */
	Eigen::RowVectorXd X_pinv = X.transpose() / (X.transpose() * X);
	Eigen::MatrixXd beta = X_pinv * Y;
	Eigen::MatrixXd resid = Y - X * beta.transpose();
	Eigen::VectorXd dof;
	if (dofPtr) {
		dof = *dofPtr;
	} else {
		dof = Eigen::VectorXd::Constant(Y.cols(), transform ? transform ->rows() : X.size());
	}

	Eigen::MatrixXd meat;
	Eigen::MatrixXd bread;
	if (!transform) {
		double XTX = X.dot(X);
		meat = resid.array().square().matrix().colwise().sum();
		bread = Eigen::MatrixXd::Constant(1, Y.cols(), 1.0 / std::max(XTX, std::numeric_limits<double>::epsilon()));

	} else {
		Eigen::VectorXd X_transpose_U = X.transpose() * (*transform).transpose();
		double X_transpose_U_U_X = X_transpose_transform.dot(*transform * X);
		meat = (transform->operator*(resid)).array().square().matrix().colwise().sum();
		bread = Eigen::MatrixXd::Constant(1, Y.cols(), 1.0 / std::max(X_transpose_U_U_X, std::numeric_limits<double>::epsilon()));
	}

	for (int i = 0; i < meat.size(); ++i) {
		meat(i) /= std::max(dof(i) - 1.0, 1.0);  // Adjusted for simple linear regression
	}

	Eigen::MatrixXd hc0_se = (bread.replicate(1, Y.cols()) * meat.array()).sqrt().matrix();
	Eigen::MatrixXd result = betas.array() / hc0_se.array();

	// Replace divisions by zero with zeros
	for (int i = 0; i < result.rows(); ++i) {
		for (int j = 0; j < result.cols(); ++j) {
			if (hc0_se(i, j) == 0) {
				result(i, j) = 0;
			}
		}
	}

	return result;
}
