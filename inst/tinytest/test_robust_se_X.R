set.seed(42)

if (requireNamespace(c("lmtest", "sandwich"), quietly=TRUE)) {
  mat_nr <- c(100, 250, 500, 1000, 2000, 5000)
  mat_nc <- sample(seq(2, length(LETTERS)), length(mat_nr))
  apply(cbind(mat_nr, mat_nc), 1, function(shape) {
    mat <- matrix(rnorm(prod(shape)), ncol=shape[["mat_nc"]], nrow=shape[["mat_nr"]])
    colnames(mat) <- LETTERS[seq_len(ncol(mat))]
    i <- sample(colnames(mat), 1)
    actual <- sapply(setdiff(colnames(mat), i), function(cn) {
      model <- lm(y~0+x, data.frame(y=mat[,cn], x=mat[,i]))
      se <- lmtest::coeftest(model, vcov.=sandwich::vcovHC(model, "HC0"))
      return(se["x", "t value"])
    })
    predicted <- scdemon::robust_se_X(i, mat)
    expect_equal(actual, predicted[names(actual)])
  })
}
