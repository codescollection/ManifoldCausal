#' Procrustes alignment and distance matrix
#'
#' @param x A k x 2 x n array of landmark configurations.
#' @return A list with components: `gpa` (procGPA result), `aligned` (aligned array),
#'   `scores` (tangent space scores), `dist` (n x n full Procrustes distance matrix).
#' @export
shape_alignment <- function(x) {
  gpa <- shapes::procGPA(x, scale = TRUE, reflect = FALSE)
  n <- dim(x)[3]
  dist_mat <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      dist_mat[i,j] <- shapes::procdist(x[,,i], x[,,j], type = "full")
      dist_mat[j,i] <- dist_mat[i,j]
    }
  }
  list(gpa = gpa, aligned = gpa$rotated, scores = gpa$scores, dist = dist_mat)
}