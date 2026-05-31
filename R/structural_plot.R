#' Plot structural function along first PC
#'
#' @param gpa A procGPA object.
#' @param gamma Structural coefficient vector (length n).
#' @param Y Outcome vector.
#' @param shape_array Original k x 2 x n array.
#' @param sigma_D Gaussian kernel bandwidth.
#' @param cities Optional city names.
#' @return A ggplot object.
#' @export
plot_structural <- function(gpa, gamma, Y, shape_array, sigma_D, cities = NULL) {
  k <- nrow(gpa$mshape)
  mshape_vec <- as.vector(gpa$mshape)
  pc1_vec <- gpa$pcar[,1]
  s_range <- range(gpa$scores[,1])
  s_seq <- seq(s_range[1]-0.5, s_range[2]+0.5, length.out=100)
  
  g_curve <- numeric(length(s_seq))
  for (i in seq_along(s_seq)) {
    s <- s_seq[i]
    new_vec <- mshape_vec + s * pc1_vec
    new_shape <- matrix(new_vec, nrow=k, ncol=2, byrow=FALSE)
    dists <- sapply(1:dim(shape_array)[3], function(j) shapes::procdist(new_shape, shape_array[,,j], type="full"))
    k_new <- gauss_kernel(dists, sigma_D)
    g_curve[i] <- sum(gamma * k_new)
  }
  
  df <- data.frame(PC1 = s_seq, g = g_curve)
  city_df <- data.frame(PC1 = gpa$scores[,1], g = Y)
  if (!is.null(cities)) city_df$city <- cities
  
  ggplot(df, aes(x = PC1, y = g)) +
    geom_line(color = "#0072B2", size = 1.2) +
    geom_point(data = city_df, aes(x = PC1, y = g), color = "#D55E00", size = 3) +
    {if (!is.null(cities)) geom_text(data = city_df, aes(label = city), vjust = -1, color = "#D55E00", size = 4.5)} +
    labs(x = "Shape PC1", y = "log GDP", title = "Structural Function: Effect of Road Network Shape on Economic Output") +
    theme_minimal(base_size = 14)
}