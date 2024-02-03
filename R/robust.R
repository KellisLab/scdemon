#' Run OLS
#'
#' @param X feature matrix
#' @param Y target
#' @return A matrix
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp evalCpp
ols_beta <- function(X, Y) {
    if (Matrix::isDiagonal(X)) {
        return(Matrix::Diagonal(x=1/Matrix::diag(X)) %*% Y)
    }
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
#' 
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
#' 
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

#' Calculate HC0 SE per-row
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp evalCpp
robust_se_X <- function(cname, Y) {
    stopifnot(cname %in% colnames(Y))
    setNames(r_robust_se_X(match(cname, colnames(Y)) - 1, Y), colnames(Y))
}


## #' Calculate HC0 SE per-row
## #' @export
## #' @useDynLib scdemon
## #' @importFrom Rcpp evalCpp
## robust_se_F <- function(cname, U, V, block_size=1000) {
##     stopifnot(cname %in% colnames(V))
##     setNames(r_robust_se_Xfull(match(cname, colnames(V)) - 1, U, V, block_size), colnames(V))
## }

#' calculate robust standard error t-values.
#' Pass U, V such that X=UV
#' @param U Observation decomposition
#' @param V Variable decomposition
#' @param B Covariates/batch effects. Uses just intercept if NULL
#' @param t_cutoff Cutoff to add to matrix. If default, uses nominal_p_cutoff
#' @param abs_t Whether to include items Pr>|t| or just Pr>t
#' @param nominal_p_cutoff Cutoff to include items for automatic filtering.
#' @return Sparse matrix of t-values, or absolute-value t-values if abs_t=T
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp evalCpp
robust_se_t.default <- function(U, V, B=NULL, t_cutoff=NULL,
                                nominal_p_cutoff=0.05,
                                abs_t=FALSE,
                                n_components=NULL,
                                min_norm=1e-300) {
    stopifnot(ncol(U)==nrow(V))
    if (!is.null(n_components)) {
        stopifnot(n_components > 0)
        U = U[ ,1:n_components, drop=FALSE]
        V = V[1:n_components, , drop=FALSE]
    }
    if (is.null(B)) {
        B = matrix(1, nrow=nrow(U))
    }
    stopifnot(nrow(U)==nrow(B))
    if (min_norm > 0) {
        V = V[, apply(V, 2, norm, "2") >= min_norm]
    }
    cat("Decomposing V with SVD\n")
    V.svd = svd(V)
    rownames(V.svd$v) = colnames(V)
    U = U %*% V.svd$u 
    cat("Extracting non-orthogonal residuals\n")
    lhs.qr = qr(ols_resid(X=B, Y=U, beta=ols_beta(X=B, Y=U)))
    cat("Computing new embedding\n")
    V = qr.R(lhs.qr) %*% diag(V.svd$d) %*% t(V.svd$v)
    if (is.null(t_cutoff)) {
        ## should be around 6.5 for most snRNA-seq datasets
        t_cutoff = qt(min(nominal_p_cutoff * ncol(V)**-2, 1),
                      nrow(B)-ncol(B),
                      lower.tail=FALSE)
    }
    M = r_robust_se(V, t_cutoff, abs_t);
    dimnames(M) = list(colnames(V), colnames(V));
    return(M);
}


#' Calculate robust standard error p-values.
#' Pass U, V such that X=UV
#' @param U Observation decomposition
#' @param V Variable decomposition
#' @param B Covariates/batch effects. Uses just intercept if NULL
#' @param nnz Number of nonzero entries per column of V, used for dof calculation.
#' @param abs_t Whether to include items Pr>|t| or just Pr>t
#' @param nominal_p_cutoff Cutoff to include items for automatic filtering.
#' @return Sparse matrix of p-values
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp evalCpp
robust_se_p.default <- function(U, V, B=NULL, nnz=NULL,
                                nominal_p_cutoff=0.05,
                                abs_t=FALSE,
                                n_components=NULL,
                                min_norm=1e-20) {
    stopifnot(ncol(U)==nrow(V))
    if (!is.null(n_components)) {
        stopifnot(n_components > 0)
        U = U[ ,1:n_components, drop=FALSE]
        V = V[1:n_components, , drop=FALSE]
    }
    if (is.null(B)) {
        B = matrix(1, nrow=nrow(U))
    }
    stopifnot(nrow(U)==nrow(B))
    if (min_norm > 0) {
        V.keep = apply(V, 2, norm, "2") >= min_norm
        V = V[, V.keep]
        if (!is.null(nnz)) {
            stopifnot(length(nnz)==length(V.keep))
            nnz = nnz[V.keep]
        }
    }
    cat("Decomposing V with SVD\n")
    V.svd = svd(V)
    rownames(V.svd$v) = colnames(V)
    U = U %*% V.svd$u 
    cat("Extracting non-orthogonal residuals\n")
    lhs.qr = qr(ols_resid(X=B, Y=U, beta=ols_beta(X=B, Y=U)))
    cat("Computing new embedding\n")
    V = qr.R(lhs.qr) %*% diag(V.svd$d) %*% t(V.svd$v)
    if (is.null(nnz)) {
        ## By default, just use all , but filter as norm approaches zero
        nnz = nrow(U) * (apply(V, 2, norm, "2") > 1e-10)
    }
    dof = pmax(nnz - ncol(B), 1)
    M = r_robust_se_p(V, dof, nominal_p_cutoff, abs_t);
    dimnames(M) = list(colnames(V), colnames(V));
    return(M);
}
