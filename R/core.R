
#' Checking PCA variance
#' @param Vh PCA loadings (PCs by features)
#' @param std Per-feature standard deviations at which the pre-PCA data matrix was used. If scaled, NULL or rep(1, ncol(Vh)) will work.
#' @param var_explained Variance explained
#' @param variance_ratio Variance ratio
#' @param tol Tolerance
#' @export
check_feature_variance <- function(Vh, std=NULL, var_explained=NULL, variance_ratio=NULL, tol=0.001) {
  stopifnot(!is.null(colnames(Vh)))
  stopifnot(ncol(Vh) == length(std))
  if (is.null(variance_ratio) || is.null(var_explained)) {
    warning("Unable to check if PCA is scaled or not, assuming correct standard deviation of features is passed...")
    stopifnot(!is.null(std))
    return(std**2)
  } 
  stopifnot(sum(variance_ratio)<=1.)
  num_nonzero_pc <- sum(apply(Vh, 2, norm, "2") > 0.)
  total_var <- sum(var_explained)/sum(variance_ratio)
  if (abs(total_var - num_nonzero_pc) < 0.001) {
    return(rep(1, ncol(Vh)))
  }
  if (is.null(std)) {
    stop("PCA is not standardized, but standard deviation used for this standardization is not passed")
  }
  if (abs(total_var - sum(std**2)) > tol) {
    warning(paste0("Standard deviation provided (", sum(std**2), ") does not match total variance (", total_var, "), continuing anyway"))
  }
  return(std**2)
}
#' SCDemon default
#' 
#' @param UD PCA embeddings (cell by PCs)
#' @param Vh PCA loadings (PC by gene)
#' @param std Per-gene standard deviation.
#' @param nominal_pvalue_cutoff Nominal Pvalue cutoff before Bonferroni
#' @param symmetric Whether to use the symmetric Wald (T stat) or asymmetric Wald (F stat)
#' @param t_cutoff Optional parameter to hard-set the T cutoff
#' @param f_cutoff Optional parameter to hard-set the F cutoff
#' @param var_explained Explained variance
#' @param variance_ratio asdf
#' @param batch_size Batching parameter to save memory
#' @return A sparse matrix containing entries with the (signed) relevant statistic if significant.
#' @export
scdemon.default <- function(UD, Vh, std, nominal_pvalue_cutoff=0.05,
                            symmetric=TRUE, t_cutoff=NA, f_cutoff=NA,
                            var_explained=NULL, variance_ratio=NULL,
                            batch_size=500L) {
  stopifnot(!is.null(colnames(Vh)))
  stopifnot(ncol(UD) == nrow(Vh))
  stopifnot(ncol(Vh) == length(std))
  stopifnot(batch_size > 0)
  stopifnot(is.integer(batch_size))
  feature_var <- check_feature_variance(Vh=Vh, std=std, var_explained=var_explained, variance_ratio=variance_ratio)
  S <- apply(UD, 2, norm, "2")
  feature_pc_var <- apply(diag(S) %*% Vh, 2, norm, "2")**2 / (nrow(UD) - 1)
  pc_vcov <- apply(diag(1/S) %*% Vh, 2, norm, "2")**2 / sqrt(sum(S**-2))
  me_var <- pmax(1e-50, feature_var - feature_pc_var)
  dof <- length(S)
  if (!all(feature_var==1.)) {
    ### Run ebayes but only on features that are captured by PCA
    eb <- ebayes(me_var[feature_pc_var > 0], d_g=dof)
    me_var[feature_pc_var > 0] <- eb
    dof <- attr(eb, "dof")
  }
  if (is.na(f_cutoff)) {
    f_cutoff <- qf(nominal_pvalue_cutoff * ncol(Vh)**-2, 2, 2*dof, lower.tail=FALSE)
  }
  if (is.na(t_cutoff)) {
    t_cutoff <- qt(nominal_pvalue_cutoff * ncol(Vh)**-2, dof, lower.tail=FALSE)
  }
  if (symmetric) {
    return(.scdemon_f(Vh=Vh, pc_vcov=pc_vcov, resid_var=me_var,
                      f_cutoff=f_cutoff, batch_size=batch_size))
  } else {
    return(.scdemon_t(Vh=Vh, pc_vcov=pc_vcov, resid_var=me_var,
                      t_cutoff=t_cutoff, batch_size=batch_size))
  }
}


#' @import Matrix
#' @importFrom progressr with_progress progressor
.scdemon_t <- function(Vh, pc_vcov, resid_var, t_cutoff, batch_size=NA) {
  if (is.nan(batch_size)) {
    batch_size <- ncol(Vh)
  }
  stopifnot(!is.null(colnames(Vh)))
  stopifnot(is.integer(batch_size))
  stopifnot(!is.nan(t_cutoff))
  stopifnot(ncol(Vh)==length(pc_vcov))
  stopifnot(length(pc_vcov)==length(resid_var))
  stopifnot(all(pc_vcov > 0))
  stopifnot(all(resid_var > 0))
  stopifnot(t_cutoff > 0)
  with_progress({
    batches <- as.integer(seq(1L, ncol(Vh), by=batch_size))
    p <- progressor(along=batches)
    do.call(cbind, lapply(batches, function(left) {
      right <- min(left + batch_size - 1, ncol(Vh))
      beta <- t(Vh) %*% Vh[,left:right]
      tval <- beta / sqrt(outer(pc_vcov, resid_var[left:right]))
      p()
      Matrix(tval * (abs(tval) >= t_cutoff),
             dimnames=list(colnames(Vh), colnames(Vh)[left:right]))
    }))
  })
}

#' @import Matrix
#' @importFrom progressr with_progress progressor
.scdemon_f <- function(Vh, pc_vcov, resid_var, f_cutoff, batch_size=NA) {
  if (is.nan(batch_size)) {
    batch_size <- ncol(Vh)
  }
  stopifnot(!is.null(colnames(Vh)))
  stopifnot(is.integer(batch_size))
  stopifnot(!is.nan(f_cutoff))
  stopifnot(ncol(Vh)==length(pc_vcov))
  stopifnot(length(pc_vcov)==length(resid_var))
  stopifnot(all(pc_vcov > 0))
  stopifnot(all(resid_var > 0))
  stopifnot(f_cutoff > 0)
  with_progress({
    batches <- as.integer(seq(1L, ncol(Vh), by=batch_size))
    p <- progressor(along=batches)
    do.call(cbind, lapply(batches, function(left) {
      right <- min(left + batch_size - 1, ncol(Vh))
      beta <- t(Vh) %*% Vh[,left:right]
      t_left <- beta / sqrt(outer(pc_vcov, resid_var[left:right]))
      t_right <- beta / sqrt(outer(resid_var, pc_vcov[left:right]))
      fstat <- abs(t_left**2 + t_right**2) / 2
      p()
      Matrix(sign(beta) * sqrt(fstat) * (fstat >= f_cutoff),
             dimnames=list(colnames(Vh), colnames(Vh)[left:right]))
    }))
  })
}
