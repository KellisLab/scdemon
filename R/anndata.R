
scdemon.AbstractAnnData <- function(adata,
                                    obsm="X_pca", varm="PCs",
                                    scaled=FALSE, nominal_pvalue_cutoff=0.05,
                                    symmetric=TRUE, t_cutoff=NA, f_cutoff=NA, batch_size=500L,
                                    inplace=TRUE, key_added="scdemon") {
  UD = adata$obsm[[obsm]]
  Vh = t(adata$varm[[varm]])
  colnames(Vh) = adata$var_names
  if (scaled) {
    std=rep(1, ncol(Vh))
  } else {
    std <- adata$var$std
  }
  M <- scdemon.default(UD, Vh, std, nominal_pvalue_cutoff=nominal_pvalue_cutoff,
                       symmetric=symmetric, t_cutoff=t_cutoff, f_cutoff=f_cutoff,
                       batch_size=batch_size)
  if (inplace) {
    adata$uns[[key_added]] = list(stat=key_added, symmetric=symmetric)
    adata$varp[[key_added]] = M
  } else {
    return(M)
  }
}
