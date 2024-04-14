#' Run OLS
#'
#' @param X feature matrix
#' @param Y target
#' @return A matrix
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp evalCpp
ols_beta <- function(X, Y, lambda=0) {
  if (Matrix::isDiagonal(X)) {
    return(Matrix::Diagonal(x=1/Matrix::diag(X)) %*% Y)
  }
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  if (!is.null(nrow(Y))) {
    stopifnot(nrow(X)==nrow(Y))
  } else {
    stopifnot(nrow(X)==length(Y))
  }
  stopifnot(lambda >= 0)
  beta <- r_ols_beta(X, Y, lambda)
  dimnames(beta) <- list(colnames(X), colnames(Y))
  return(beta)
}

#' Calculate OLS residuals
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp evalCpp
ols_resid <- function(X, Y, beta) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  beta <- as.matrix(beta)
  stopifnot(nrow(X)==nrow(Y))
  stopifnot(ncol(X)==nrow(beta))
  stopifnot(ncol(beta)==ncol(Y))
  res <- r_ols_resid(X, Y, beta)
  dimnames(res) <- dimnames(Y)
  return(res)
}

#' Calculate HC0 SE per-row
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp evalCpp
robust_se_X <- function(cname, Y, lambda=1e-10) {
  stopifnot(cname %in% colnames(Y))
  setNames(r_robust_se_X(Y[,match(cname, colnames(Y))], Y, lambda), colnames(Y))
}

#' Calculate regularized SE with per-Y lambda
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp evalCpp
robust_se_L <- function(cname, Y, lambda) {
    stopifnot(cname %in% colnames(Y))
    stopifnot(length(lambda)==ncol(Y))
    setNames(r_robust_se_L(Y[,match(cname, colnames(Y))], Y, lambda), colnames(Y))
}

### TODO test Y=U when U is diagonal for ols_resid. Only missing part for use=X
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp evalCpp
.robust_prepare <- function(U, V, B=NULL, n_components=NULL, min_norm=1e-5, return_U=FALSE, lambda=0, power=1) {
  stopifnot(ncol(U)==nrow(V))
  if (!is.null(n_components)) {
    stopifnot(n_components > 0)
    U <- U[ ,seq_len(n_components), drop=FALSE]
    V <- V[seq_len(n_components), , drop=FALSE]
  }
  if (is.null(B)) {
    B <- matrix(1, nrow=nrow(U))
  }
  stopifnot(nrow(U)==nrow(B))
  rownames(B) <- rownames(U)
  if (min_norm > 0) {
    cat("Filtering norm\n")
    V <- V[, apply(V, 2, norm, "2") >= min_norm]
  }
  xform = ols_beta(X=U, Y=B, lambda=lambda) %*% ols_beta(X=B, Y=U, lambda=lambda)
  cat("Decomposing V with SVD\n")
  V_svd <- svd(V - xform %*% V)
  rownames(V_svd$v) <- colnames(V)
  ## group orthogonal items
  ## TODO sample sqrt(U)? Experiment with U sampling size

  cat("Extracting non-orthogonal residuals\n")
  QR <- qr(U %*% V_svd$u)
  cat("Computing new embedding\n")
  RS_svd <- svd(qr.R(QR) %*% diag(V_svd$d))
  V <- diag(RS_svd$d**power) %*% t(RS_svd$v) %*% t(V_svd$v)
  attr(V, "dof") <- nrow(B) - ncol(B)
  if (return_U) {
    U <- qr.Q(QR) %*% RS_svd$u %*% diag(RS_svd$d ** (1-power))
    rownames(U) <- rownames(B)
    attr(V, "U") <- U
  }
  return(V)
}
#' calculate robust standard error t-values.
#' Pass U, V such that X=UV
#' @param V De-biased 
#' @param V Variable decomposition
#' @param B Covariates/batch effects. Uses just intercept if NULL
#' @param t_cutoff Cutoff to add to matrix. If default, uses nominal_p_cutoff
#' @param abs_t Whether to include items Pr>|t| or just Pr>t
#' @param nominal_p_cutoff Cutoff to include items for automatic filtering.
#' @return Sparse matrix of t-values, or absolute-value t-values if abs_t=T
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp evalCpp
robust_se_t.default <- function(Us, V,
                                nominal_p_cutoff=0.05,
                                lambda=100,
                                abs_t=FALSE, dof=NULL,
                                t_cutoff=NULL) {
  require(Matrix)
  if (is.null(t_cutoff)) {
    ## should be around 6.5 for most snRNA-seq datasets
    t_cutoff <- qt(min(nominal_p_cutoff / (ncol(V) * ncol(V)), 1),
                   max(nrow(Us) - ncol(V), 4),
                   lower.tail=FALSE)
  }
  M <- r_robust_se(Us, V, lambda, t_cutoff, abs_t)
  dimnames(M) <- list(colnames(V), colnames(V))
  return(Matrix::t(Matrix::drop0(M)))
}

