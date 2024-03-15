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
.robust_prepare <- function(U, V, B=NULL, n_components=NULL, min_norm=1e-5, return_U=FALSE) {
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
  cat("Decomposing V with SVD\n")
  V_svd <- svd(V)
  rownames(V_svd$v) <- colnames(V)
  ## group orthogonal items
  ## TODO sample sqrt(U)? Experiment with U sampling size
  U <- U %*% V_svd$u 
  cat("Extracting non-orthogonal residuals\n")
  lhs_qr <- qr(ols_resid(X=B, Y=U, beta=ols_beta(X=B, Y=U)))
  cat("Computing new embedding\n")
  RS_svd <- svd(qr.R(lhs_qr) %*% diag(V_svd$d))
  V <- diag(RS_svd$d) %*% t(RS_svd$v) %*% t(V_svd$v)
  attr(V, "dof") <- nrow(B) - ncol(B)
  if (return_U) {
    U <- qr.Q(lhs_qr) %*% RS_svd$u
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
robust_se_t.default <- function(V1, V2,
                                nominal_p_cutoff=0.05,
                                lambda=1e-10,
                                abs_t=FALSE, dof=NULL,
                                t_cutoff=NULL) {
  require(Matrix)
  if (is.null(V2)) V2 <- V1 
  if (is.null(t_cutoff)) {
    if (is.null(dof)) {
      stopifnot(is.integer(attr(V1, "dof")))
      stopifnot(is.integer(attr(V2, "dof")))
      dof <- mean(c(attr(V1, "dof"), attr(V2, "dof")))
    }
    ## should be around 6.5 for most snRNA-seq datasets
    t_cutoff <- qt(min(nominal_p_cutoff / (ncol(V1) * ncol(V2)), 1),
                   dof,
                   lower.tail=FALSE)
  }
  comm <- intersect(colnames(V2), colnames(V1))
  M <- r_robust_se(V1, V2, lambda, t_cutoff, abs_t)
  dimnames(M) <- list(colnames(V2), colnames(V1))
  if (!is.null(dof)) attr(M, "dof") <- dof
  if (length(comm) > 0) { 
    M[cbind(match(colnames(V2), comm),
            match(colnames(V1), comm))] <- 0
  }
  return(Matrix::t(Matrix::drop0(M)))
}

