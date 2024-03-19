// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// r_ols_beta
Eigen::MatrixXd r_ols_beta(const Eigen::Map<Eigen::MatrixXd>& X, const Eigen::Map<Eigen::MatrixXd>& Y, double lambda);
RcppExport SEXP _scdemon_r_ols_beta(SEXP XSEXP, SEXP YSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(r_ols_beta(X, Y, lambda));
    return rcpp_result_gen;
END_RCPP
}
// r_ols_resid
Eigen::MatrixXd r_ols_resid(const Eigen::Map<Eigen::MatrixXd>& X, const Eigen::Map<Eigen::MatrixXd>& Y, const Eigen::Map<Eigen::MatrixXd>& beta);
RcppExport SEXP _scdemon_r_ols_resid(SEXP XSEXP, SEXP YSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(r_ols_resid(X, Y, beta));
    return rcpp_result_gen;
END_RCPP
}
// r_robust_se_X
Eigen::ArrayXd r_robust_se_X(const Eigen::Map<Eigen::MatrixXd>& X, const Eigen::Map<Eigen::MatrixXd>& Y, double lambda);
RcppExport SEXP _scdemon_r_robust_se_X(SEXP XSEXP, SEXP YSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(r_robust_se_X(X, Y, lambda));
    return rcpp_result_gen;
END_RCPP
}
// r_robust_se_Y
Eigen::ArrayXd r_robust_se_Y(const Eigen::Map<Eigen::MatrixXd>& X, const Eigen::Map<Eigen::MatrixXd>& Y, double lambda);
RcppExport SEXP _scdemon_r_robust_se_Y(SEXP XSEXP, SEXP YSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(r_robust_se_Y(X, Y, lambda));
    return rcpp_result_gen;
END_RCPP
}
// r_ols_beta_L
Eigen::MatrixXd r_ols_beta_L(const Eigen::Map<Eigen::VectorXd>& X, const Eigen::Map<Eigen::MatrixXd>& Y, const Eigen::Map<Eigen::ArrayXd>& lambda);
RcppExport SEXP _scdemon_r_ols_beta_L(SEXP XSEXP, SEXP YSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXd>& >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(r_ols_beta_L(X, Y, lambda));
    return rcpp_result_gen;
END_RCPP
}
// r_robust_se_L
Eigen::ArrayXd r_robust_se_L(const Eigen::Map<Eigen::MatrixXd>& X, const Eigen::Map<Eigen::MatrixXd>& Y, const Eigen::Map<Eigen::ArrayXd>& lambda);
RcppExport SEXP _scdemon_r_robust_se_L(SEXP XSEXP, SEXP YSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXd>& >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(r_robust_se_L(X, Y, lambda));
    return rcpp_result_gen;
END_RCPP
}
// r_robust_se
Eigen::SparseMatrix<double> r_robust_se(const Eigen::Map<Eigen::MatrixXd>& X, const Eigen::Map<Eigen::MatrixXd>& Y, double lambda, double t_cutoff, bool abs_cutoff);
RcppExport SEXP _scdemon_r_robust_se(SEXP XSEXP, SEXP YSEXP, SEXP lambdaSEXP, SEXP t_cutoffSEXP, SEXP abs_cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type t_cutoff(t_cutoffSEXP);
    Rcpp::traits::input_parameter< bool >::type abs_cutoff(abs_cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(r_robust_se(X, Y, lambda, t_cutoff, abs_cutoff));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scdemon_r_ols_beta", (DL_FUNC) &_scdemon_r_ols_beta, 3},
    {"_scdemon_r_ols_resid", (DL_FUNC) &_scdemon_r_ols_resid, 3},
    {"_scdemon_r_robust_se_X", (DL_FUNC) &_scdemon_r_robust_se_X, 3},
    {"_scdemon_r_robust_se_Y", (DL_FUNC) &_scdemon_r_robust_se_Y, 3},
    {"_scdemon_r_ols_beta_L", (DL_FUNC) &_scdemon_r_ols_beta_L, 3},
    {"_scdemon_r_robust_se_L", (DL_FUNC) &_scdemon_r_robust_se_L, 3},
    {"_scdemon_r_robust_se", (DL_FUNC) &_scdemon_r_robust_se, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_scdemon(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
