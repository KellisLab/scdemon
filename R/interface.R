#' Run OLS
#'
#' @param a vector
#' @param b matrix
#' @return A matrix
#' @export
#' @useDynLib scdemon
#' @importFrom Rcpp sourceCpp
ols <- function(a, b) {
    return(r_ols(a, b))
}
