#' Confidence ellipse and outliers from PC1/PC2
#'
#' For computing a confidence ellipse for spectra pc1/pc2
#' can use with either raw or processed but recommend processed
#' since most of the variance occurs in the intensities
#' @param spectra Spectra data that has been interpolated already
#' @param confidence Chosen confidence level, recommend 95%
#' @return List of 3 elements: Matrix with PC1/PC2 of your data,
#' Confidence ellipse matrix, the indices of any outliers
#' @export
#' @import pracma
conf_ellipse_outliers <- function(spectra, confidence){
  # Compute PCA on spectra
  spec_pca <- prcomp(spectra, center = TRUE, scale = TRUE)

  # Pick off only PC1 and PC2
  mat <- as.matrix(spec_pca)

  # Mean matrix
  mean_mat <- apply(mat, 2, mean)

  # Compute covariance matrix
  cov_mat <- cov(mat)

  # Eigenvectors and Eigenvalues for the covariance matrix
  evals <- eigen(cov_mat)$values
  evecs <- eigen(cov_mat)$vectors

  # Angles of a circle
  a <- seq(0, 2*pi, len=1000)

  # Get critical value
  c2 <- qchisq(as.numeric(confidence), 2) # uses chisquare distribution set to X confidence
  c <- sqrt(c2)

  # Get the distances
  xT <- c * sqrt(evals[1]) * cos(a)
  yT <- c * sqrt(evals[2]) * sin(a)

  M <- cbind(xT, yT)

  # Covert the coordinates
  transM <- evecs %*% t(M)
  transM <- t(transM)

  # Add back on the mean to centre
  output <- transM + mean_mat

  # Store indices of those which lie outside of the X confidence ellipse
  outliers <- which(inpolygon(x = mat[,1], mat[,2], transM[,1], transM[,2], boundary = FALSE) == FALSE)

  return(list(mat, output, outliers))
}
