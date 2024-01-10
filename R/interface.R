#' Run OLS
#'
#' @param X feature matrix
#' @param Y target
#' @return A matrix
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp evalCpp
ols_beta <- function(X, Y) {
    X = as.matrix(X)
    if (!is.null(nrow(Y))) {
        stopifnot(nrow(X)==nrow(Y))
    } else {
        stopifnot(nrow(X)==length(Y))
    }
    return(r_ols_beta(X, Y))
}

#' Calculate OLS residuals
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp evalCpp
ols_resid <- function(X, Y, beta) {
    X = as.matrix(X)
    Y = as.matrix(Y)
    beta = as.matrix(beta)
    stopifnot(nrow(X)==nrow(Y))
    stopifnot(ncol(X)==nrow(beta))
    stopifnot(ncol(beta)==ncol(Y))
    return(r_ols_resid(X, Y, beta))
}

#' Adjust residuals using FW partialling out
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp evalCpp
fw_meat <- function(res, U=NULL, B=NULL, BPU=NULL) {
    if (is.null(U)) {
        return(res)
    }
    if (!is.matrix(res)) {
        res = as.matrix(res)
    }
    stopifnot(is.matrix(U))
    stopifnot(ncol(U)==nrow(res))
    if (is.null(B)) {
        B = matrix(1, nrow=nrow(U), ncol=1)
    }
    stopifnot(is.matrix(B))
    stopifnot(nrow(B)==nrow(U))
    if (is.null(BPU)) {
        BPU = MASS::ginv(B) %*% U
    }
    stopifnot(is.matrix(BPU))
    stopifnot(nrow(BPU)==ncol(B))
    stopifnot(ncol(BPU)==ncol(U))
    return(r_fw_meat(res, U, B, BPU))
}


#' Adjust X using FW partialling out
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp evalCpp
fw_bread <- function(X, U=NULL, B=NULL, BPU=NULL) {
    if (!is.matrix(X)) {
        X = as.matrix(X)
    }
    if (is.null(U)) {
        U = diag(nrow(X))
    }
    stopifnot(is.matrix(U))
    stopifnot(ncol(U)==nrow(X))
    if (is.null(B)) {
        B = matrix(1, nrow=nrow(U), ncol=1)
    }
    stopifnot(is.matrix(B))
    stopifnot(nrow(B)==nrow(U))
    if (is.null(BPU)) {
        BPU = MASS::ginv(B) %*% U
    }
    stopifnot(is.matrix(BPU))
    stopifnot(nrow(BPU)==ncol(B))
    stopifnot(ncol(BPU)==ncol(U))
    return(r_fw_bread(X, U, B, BPU))
}
#' Calculate HC0 SE using U as a transform, and B as a covariate matrix.
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp evalCpp
hc0_se_Xvec <- function(X, Y, U=NULL, B=NULL, BPU=NULL) {
    stopifnot(nrow(X)==nrow(Y))
    if (is.null(U)) {
        U = diag(nrow(Y))
    }
    stopifnot(ncol(U)==nrow(Y))
    if (is.null(B)) {
        B = matrix(1, nrow=nrow(U), ncol=1)
    }
    stopifnot(is.matrix(B))
    stopifnot(nrow(B)==nrow(U))
    if (is.null(BPU)) {
        BPU = MASS::ginv(B) %*% U
    }
    stopifnot(is.matrix(BPU))
    stopifnot(nrow(BPU)==ncol(B))
    stopifnot(ncol(BPU)==ncol(U))
    return(r_hc0_se_Xvec(X, Y, U, B, BPU))
}
