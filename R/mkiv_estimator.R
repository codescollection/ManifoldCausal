#' Two-stage manifold kernel IV estimator
#'
#' @param Y Numeric vector of outcomes (length n).
#' @param K_D Kernel matrix (n x n) of shape treatment.
#' @param K_Z Kernel matrix (n x n) of instrument.
#' @param lambda Regularization parameter for the second stage.
#' @param rho1 Regularization for the first stage.
#' @return A list with `fitted` (fitted values), `alpha` (second-stage coefficients),
#'   `gamma` (structural function coefficients on original kernel).
#' @export
mkiv_fit <- function(Y, K_D, K_Z, lambda = 0.1, rho1 = 0.1) {
  n <- length(Y)
  # Stage 1: predict shape kernel from instrument kernel
  A <- chol2inv(chol(K_Z + rho1 * diag(n))) %*% K_D
  K_D_hat <- K_Z %*% A
  # Stage 2: regress Y on predicted kernel
  alpha <- chol2inv(chol(t(K_D_hat) %*% K_D_hat + lambda * diag(n))) %*% t(K_D_hat) %*% Y
  fitted <- K_D_hat %*% alpha
  # Structural function coefficients gamma = M %*% alpha, where M = K_Z (K_Z + rho I)^{-1}
  M <- K_Z %*% chol2inv(chol(K_Z + rho1 * diag(n)))
  gamma <- M %*% alpha
  list(fitted = fitted, alpha = alpha, gamma = gamma)
}