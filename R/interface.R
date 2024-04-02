
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

#' 
#' @export
#' @method robust_se_t Seurat
robust_se_t.Seurat <- function(obj, covariates=NULL,
                               reduction="pca", ### todo multiomic for multiple reductions?
                               key_added="scdemon",
                               nominal_p_cutoff=0.05,
                               t_cutoff=NULL, abs_t=FALSE, n_components=NULL,
                               min_norm=1e-5) {
  require(SeuratObject)
  U <- Embeddings(obj, reduction=reduction)
  V <- t(Loadings(obj, reduction=reduction))
  B <- .extract_covariates(covariates, df=obj[[]])
  V <- .robust_prepare(U=U, V=V, B=B, n_components=n_components, min_norm=min_norm, return_U=FALSE)
  S <- as.Graph(robust_se_t.default(V, V, t_cutoff=t_cutoff,
                                    abs_t=abs_t, nominal_p_cutoff=nominal_p_cutoff))
  slot(object = S, name = "assay.used") <- DefaultAssay(object=obj[[reduction]])
  ## TODO add to object
  return(S)
}

#' @export
robust_se_t.AbstractAnnData <- function(obj, covariates=NULL,
                                        method="pca",
                                        key_added="scdemon",
                                        nominal_p_cutoff=0.05,
                                        t_cutoff=NULL, abs_t=FALSE, n_components=NULL,
                                        min_norm=1e-5) {
  if (length(method) > 1) {
    U <- obj$obsm[[method[[1]] ]]
    V <- t(obj$varm[[method[[2]] ]])
  } else if (method == "pca") {
    U <- obj$obsm$X_pca
    V <- t(obj$varm$PCs)
  } else if (method == "lsi") {
    U <- obj$obsm$X_lsi
    V <- t(obj$varm$LSI)
  } else if (method == "X") {
    V <- obj$X
    U <- Matrix::Diagonal(n=nrow(V))
  } else if (method %in% names(obj$layers)) {
    V <- obj$layers[[method]]
    U <- Matrix::Diagonal(n=nrow(V))
  } else {
    stop(paste0("Unknown method ", paste0(method, collapse=" ")))
  }
  colnames(V) <- obj$var_names
  rownames(U) <- obj$obs_names
  B <- .extract_covariates(covariates, obj$obs)
  V <- .robust_prepare(U=U, V=V, B=B, n_components=n_components, min_norm=min_norm, return_U=FALSE)
  S <- robust_se_t.default(V, V, lambda=1e-10, t_cutoff=t_cutoff,
                           abs_t=abs_t, nominal_p_cutoff=nominal_p_cutoff)
  dimnames(S) = list(colnames(V), colnames(V))
  if (nrow(S) != length(obj$var_names)) {
    D <- Matrix::sparseMatrix(i=match(rownames(S), obj$var_names),
                              j=seq_len(nrow(S)),
                              dims=c(length(obj$var_names),
                                     nrow(S)),
                              dimnames=list(obj$var_names, NULL))
    S <- D %*% S %*% Matrix::t(D)
  }
  obj$varp[[key_added]] <- S
  return(obj)
}

robust_se_t.MultiAssayExperiment <- function(obj, covariates=NULL) {

}
