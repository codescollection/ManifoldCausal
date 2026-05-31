#' Leave-one-out cross validation for MKIV
#'
#' @param Y Outcome vector.
#' @param K_D Shape kernel matrix.
#' @param K_Z Instrument kernel matrix.
#' @param lambda_grid Grid of lambda values to search.
#' @param rho1 First-stage regularization.
#' @return Optimal lambda.
#' @export
cv_mkiv <- function(Y, K_D, K_Z, lambda_grid = 10^seq(-3, 1, length.out=20), rho1 = 0.1) {
  n <- length(Y)
  errors <- sapply(lambda_grid, function(lam) {
    err <- 0
    for (i in 1:n) {
      train <- setdiff(1:n, i)
      K_D_train <- K_D[train, train]
      K_Z_train <- K_Z[train, train]
      Y_train <- Y[train]
      fit <- tryCatch(mkiv_fit(Y_train, K_D_train, K_Z_train, lambda = lam, rho1 = rho1),
                      error = function(e) list(fitted = rep(NA, n-1)))
      if (any(is.na(fit$fitted))) return(Inf)
      k_D_test <- K_D[i, train]
      k_Z_test <- K_Z[i, train]
      M_train <- K_Z_train %*% chol2inv(chol(K_Z_train + rho1 * diag(n-1)))
      gamma_train <- M_train %*% fit$alpha
      pred <- sum(gamma_train * k_D_test)
      err <- err + (Y[i] - pred)^2
    }
    sqrt(err/n)
  })
  lambda_grid[which.min(errors)]
}