
#' @export
robust_se <- function(obj, ...) {
    UseMethod(generic="robust_se", object=obj)
}

.extract_covariates <- function(covariates, df) {
    if (is.null(covariates)) {
        return(model.matrix(~1, data=df))
    } else if (is.formula(covariates)) {
        return(model.matrix(covariates, data=df))
    } else if (is.vector(covariates)) {
        return(model.matrix(as.formula(paste0("~", paste0(covariates, collapse="+"))), df))
    } else if (is.matrix(covariates)) {
        return(covariates)
    } else {
        cat("Covariates are unknown form. Returning empty model.\n")
        return(model.matrix(~1, data=df))
    }
}

#' @export
robust_se.AbstractAnnData <- function(obj, method="pca",
                                      key_added="scdemon",
                                      covariates=NULL,
                                      nominal_p_cutoff=0.05,
                                      t_cutoff=NULL, abs_t=FALSE) {
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
    S = robust_se.default(U=U, V=V, B=B, t_cutoff=t_cutoff,
                          abs_t=abs_t, nominal_p_cutoff=nominal_p_cutoff)
    adata$varp[[key_added]] = S
    return(adata)
}

