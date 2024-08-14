

#' @importFrom S4Vectors SelfHits mcols
#' @export
scdemon.SingleCellExperiment <- function(sce, std, reduction="PCA",
                                         nominal_pvalue_cutoff=0.05,
                                         inplace=TRUE, key_added="scdemon",
                                         t_cutoff=NA, f_cutoff=NA, symmetric=TRUE, batch_size=500L) {
  require(SingleCellExperiment)
  UD <- reducedDim(sce, reduction)
  Vh <- t(attr(UD, "rotation"))
  variance_ratio <- attr(UD, "percentVar")/100
  var_explained <- attr(UD, "varExplained")
  W <- scdemon.default(UD=UD, Vh=Vh, std=rowData(sce)[[std]],
                       nominal_pvalue_cutoff=nominal_pvalue_cutoff,
                       t_cutoff=t_cutoff, f_cutoff=f_cutoff,
                       symmetric=symmetric, batch_size=batch_size,
                       var_explained=var_explained,
                       variance_ratio=variance_ratio)
  if (inplace) {
    wf <- Matrix::summary(W)
    hits <- SelfHits(from=wf$i, to=wf$j, nnode=nrow(sce))
    mcols(hits)$stat <- wf$x
    rowPair(sce, key_added) <- hits
    return(sce)
  } else {
    return(W)
  }
}
