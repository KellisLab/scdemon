
#' @export
robust_se_p <- function(obj, ...) {
    UseMethod(generic="robust_se_p", object=obj)
}

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
    if (length(method) == 2) {
        U = obj$obsm[[method[[1]] ]]
        V = obj$varm[[method[[2]] ]]
        u2norm = apply(U, 2, norm, "2")
        U = U %*% diag(x=1/u2norm)
        V = diag(x=u2norm) %*% t(V)
    } else if (method == "pca") {
        s = sqrt(obj$uns$pca$variance * (obj$n_obs()-1))
        U = obj$obsm$X_pca %*% diag(x=1/s)
        V = diag(x=s) %*% t(obj$varm$PCs)
    } else if (method == "lsi") {
        s = obj$uns$lsi$stdev * sqrt(obj$n_obs()-1)
        U = obj$obsm$X_lsi %*% diag(x=1/s)
        V = diag(x=s) %*% t(obj$varm$LSI)
    } else {
        stop(paste0("Unknown method ", paste0(method, collapse=" ")))
    }
    colnames(V) = obj$var_names
    rownames(U) = obj$obs_names
    B = .extract_covariates(covariates, obj$obs)
    S = robust_se_t.default(U=U, V=V, B=B, t_cutoff=t_cutoff,
                            abs_t=abs_t, nominal_p_cutoff=nominal_p_cutoff,
                            n_components=n_components)
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
    if (length(method) == 2) {
        U = obj$obsm[[method[[1]] ]]
        V = obj$varm[[method[[2]] ]]
        u2norm = apply(U, 2, norm, "2")
        U = U %*% diag(x=1/u2norm)
        V = diag(x=u2norm) %*% t(V)
    } else if (method == "pca") {
        s = sqrt(obj$uns$pca$variance * (obj$n_obs()-1))
        U = obj$obsm$X_pca %*% diag(x=1/s)
        V = diag(x=s) %*% t(obj$varm$PCs)
    } else if (method == "lsi") {
        s = obj$uns$lsi$stdev * sqrt(obj$n_obs()-1)
        U = obj$obsm$X_lsi %*% diag(x=1/s)
        V = diag(x=s) %*% t(obj$varm$LSI)
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
    adata$varp[[key_added]] = S
    return(adata)
}
