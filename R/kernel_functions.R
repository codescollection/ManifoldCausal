#' Gaussian kernel
gauss_kernel <- function(D, sigma) {
  exp(- (D^2) / (2 * sigma^2))
}