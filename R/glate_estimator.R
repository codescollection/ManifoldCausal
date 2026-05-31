#' Estimate Geodesic Local Average Treatment Effect (G-LATE)
#'
#' @param Y Array k x 2 x n of shape outcomes.
#' @param T Binary treatment vector.
#' @param Z Binary instrument vector.
#' @return Estimated G-LATE distance.
#' @export
glate <- function(Y, T, Z) {
  # Compute weights for compliers (standard LATE formula)
  pZ1 <- mean(Z == 1)
  pZ0 <- mean(Z == 0)
  pT1_given_Z1 <- mean(T[Z == 1] == 1)
  pT1_given_Z0 <- mean(T[Z == 0] == 1)
  delta <- pT1_given_Z1 - pT1_given_Z0
  pi1 <- pZ1 / delta
  pi0 <- pZ0 / delta
  
  # Positive subpopulation: Z=1,T=1; negative: Z=0,T=1
  idx11 <- which(Z == 1 & T == 1)
  idx01 <- which(Z == 0 & T == 1)
  shapes_1 <- Y[,,c(idx11, idx01)]
  w_1 <- c(rep(pi1, length(idx11)), rep(-pi0, length(idx01)))
  
  # Weighted Fréchet mean for Y(1) compliers (simplified tangent space weighted mean)
  gpa1 <- shapes::procGPA(shapes_1, scale = TRUE, reflect = FALSE)
  w_norm1 <- w_1 / sum(w_1)
  weighted_scores1 <- colSums(gpa1$scores * w_norm1)
  mu1 <- gpa1$mshape + matrix(weighted_scores1[1:(2*nrow(gpa1$mshape)-4)] %*% t(gpa1$pcar[,1:(2*nrow(gpa1$mshape)-4)]), ncol=2) # simplified
  # Simplified approach: just use weighted average of aligned coordinates
  mu1 <- apply(gpa1$rotated[,,], c(1,2), function(x) sum(x * w_norm1))
  
  # Similarly for Y(0) compliers: Z=1,T=0 and Z=0,T=0
  idx10 <- which(Z == 1 & T == 0)
  idx00 <- which(Z == 0 & T == 0)
  shapes_0 <- Y[,,c(idx10, idx00)]
  w_0 <- c(rep(pi1, length(idx10)), rep(-pi0, length(idx00)))
  gpa0 <- shapes::procGPA(shapes_0, scale = TRUE, reflect = FALSE)
  w_norm0 <- w_0 / sum(w_0)
  mu0 <- apply(gpa0$rotated[,,], c(1,2), function(x) sum(x * w_norm0))
  
  # Geodesic distance between complier Fréchet means
  shapes::procdist(mu1, mu0, type = "full")
}