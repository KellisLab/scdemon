
#' @export
robust_se_p <- function(obj, ...) {
    UseMethod(generic="robust_se_p", object=obj)
}

#' @export
robust_se_t <- function(obj, ...) {
    UseMethod(generic="robust_se_t", object=obj)
}
.extract_covariates <- function(covariates, df) {
    if (is.null(covariates)) {
        return(model.matrix(~1, data=df))
    } else if (class(covariates)=="formula") {
        return(model.matrix(covariates, data=df))
    } else if (class(covariates)=="character") {
        return(model.matrix(as.formula(paste0("~", paste0(covariates, collapse="+"))), df))
    } else if (is.matrix(covariates)) {
        return(covariates)
    } else {
        cat("Covariates are unknown form. Returning empty model.\n")
        return(model.matrix(~1, data=df))
    }
}

#' @export
robust_se_t.AbstractAnnData <- function(obj, covariates=NULL,
                                      method="pca",
                                      key_added="scdemon",
                                      nominal_p_cutoff=0.05,
                                      t_cutoff=NULL, abs_t=FALSE, n_components=NULL) {
    if (length(method) > 1) {
        U = obj$obsm[[method[[1]] ]]
        V = obj$varm[[method[[2]] ]]
    } else if (method == "pca") {
        U = obj$obsm$X_pca
        V = t(obj$varm$PCs)
    } else if (method == "lsi") {
        U = obj$obsm$X_lsi
        V = t(obj$varm$LSI)
    } else if (method == "X") {
        V = obj$X
        U = Matrix::Diagonal(n=nrow(V))
    } else if (method %in% names(obj$layers)) {
        V = obj$layers[[method]]
        U = Matrix::Diagonal(n=nrow(V))
    } else {
        stop(paste0("Unknown method ", paste0(method, collapse=" ")))
    }
    colnames(V) = obj$var_names
    rownames(U) = obj$obs_names
    B = .extract_covariates(covariates, obj$obs)
    S = robust_se_t.default(U=U, V=V, B=B, t_cutoff=t_cutoff,
                            abs_t=abs_t, nominal_p_cutoff=nominal_p_cutoff,
                            n_components=n_components)
    if (nrow(S) != length(adata$var_names)) {
        ### Expand 
        D = Matrix::sparseMatrix(
                        i=match(rownames(S), adata$var_names),
                        j=1:nrow(S),
                        dims=c(length(adata$var_names),
                               nrow(S)),
                        dimnames=list(adata$var_names, NULL))
        S = D %*% S %*% Matrix::t(D)
    }

    adata$varp[[key_added]] = S
    return(adata)
}

#' @export
robust_se_p.AbstractAnnData <- function(obj, covariates=NULL,
                                        method="pca",
                                        key_added="scdemon",
                                        nominal_p_cutoff=0.05,
                                        nnz="n_cells_by_counts", abs_t=FALSE,
                                        n_components=NULL) {
    if (length(method) > 1) {
        U = obj$obsm[[method[[1]] ]]
        V = obj$varm[[method[[2]] ]]
    } else if (method == "pca") {
        U = obj$obsm$X_pca
        V = t(obj$varm$PCs)
    } else if (method == "lsi") {
        U = obj$obsm$X_lsi
        V = t(obj$varm$LSI)
    } else if (method == "X") {
        V = obj$X
        U = Matrix::Diagonal(n=nrow(V))
    } else if (method %in% names(obj$layers)) {
        V = obj$layers[[method]]
        U = Matrix::Diagonal(n=nrow(V))
    } else {
        stop(paste0("Unknown method ", paste0(method, collapse=" ")))
    }
    colnames(V) = obj$var_names
    rownames(U) = obj$obs_names
    if (nnz %in% names(obj$var)) {
        nnz = obj$var[[nnz]] 
    } else {
        nnz = NULL
    }
    B = .extract_covariates(covariates, obj$obs)
    S = robust_se_p.default(U=U, V=V, B=B, nnz=nnz,
                            abs_t=abs_t, nominal_p_cutoff=nominal_p_cutoff,
                            n_components=n_components)
    if (nrow(S) != length(adata$var_names)) {
        ### Expand 
        D = Matrix::sparseMatrix(
                        i=match(rownames(S), adata$var_names),
                        j=1:nrow(S),
                        dims=c(length(adata$var_names),
                               nrow(S)),
                        dimnames=list(adata$var_names, NULL))
        S = D %*% S %*% Matrix::t(D)
    }
    adata$varp[[key_added]] = S
    return(adata)
}
