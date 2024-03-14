set.seed(0)

mat_nr <- c(100, 250, 500, 1000, 2000, 5000)
mat_nc <- c(10, 1000, 400, 250, 555, 100)
mat_nm <- as.integer(ceiling(sqrt(pmin(mat_nr, mat_nc))))
mat_nb <- c(0, 1, 5, 10, 3, 50)
apply(cbind(mat_nr, mat_nc, mat_nm, mat_nb), 1, function(shape) {
    U <- matrix(rnorm(shape[["mat_nr"]] * shape[["mat_nm"]]), nrow=shape[["mat_nr"]])
    V <- matrix(rnorm(shape[["mat_nm"]] * shape[["mat_nc"]]), ncol=shape[["mat_nc"]])
    if (shape[["mat_nb"]] > 0) {
        ### Make B such that U has a substantial component made from B
        B <- U %*% matrix(rnorm(shape[["mat_nm"]] * shape[["mat_nb"]]), ncol=shape[["mat_nb"]])
        B <- B + matrix(rnorm(prod(dim(B))), ncol=ncol(B))
    } else {
        B <- matrix(1, nrow=shape[["mat_nr"]])
    }
    #U <- U + 5 *  rowMeans(B) ### to add effect
    ### Test that B is correctly partialled out
    uvb <- U %*% V - B %*% solve(t(B) %*% B) %*% t(B) %*% U %*% V
    V1 <- scdemon:::.robust_prepare(U, V, B, return_U=TRUE)
    U1 <- attr(V1, "U")
    expect_equal(U1 %*% V1, uvb)
    ## plot(atanh(cor(V1)), atanh(cor(uvb)))
### test accuracy
    ## I = sample(nrow(U), as.integer(ceiling(sqrt(nrow(U)))), replace=FALSE)
    ## V2 <- scdemon:::.robust_prepare(U[I,], V, B[I,], return_U=TRUE)
    ## U2 <- attr(V2, "U")
    ## plot(atanh(cor(V2)), atanh(cor(uvb)))
})
