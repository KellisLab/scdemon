
#' Empirical Bayes
#'
#' Empirical bayes to stabilize variance similar to Limma-trend, but use
#' median/MAD instead of mean/sd for robustness
#' 
#' @param s_g_squared Residual variance input
#' @param d_g Degrees of freedom for residual variance
#' @param d_0 Optional null dof, shouldn't be used 
#' @importFrom limma trigammaInverse
#' @return Array of smoothed variances, with a dof attribute containing smoothed degrees of freedom
#' @export
ebayes <- function(s_g_squared, d_g, d_0=NA) {
  n <- length(s_g_squared)
  e_g <- log(s_g_squared) - digamma(d_g/2) + log(d_g / 2)
  e_bar <- median(e_g)
  tgi_in <- median((e_g - e_bar)**2 * n/(n-1) - trigamma(d_g / 2))
  if (is.na(d_0)) {
    d_0 <- 2 * trigammaInverse(tgi_in)
  }
  s_0_squared <- exp(e_bar + digamma(d_0 / 2) - log(d_0 / 2))
  smoothed_var <- (d_0 * s_0_squared + d_g * s_g_squared)/(d_0 + d_g)
  attr(smoothed_var, "dof") <- d_g + d_0
  return(smoothed_var)
}
