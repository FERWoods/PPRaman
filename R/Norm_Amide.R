#' Normalise at the largest amide peak
#'
#' Normalising function using min/max method at amide peak
#' @param spectra Spectra in matrix format for normalising
#' @return Normalised spectra
#' @export
norm_a <- function(spectra){
  norm_spectra <- matrix(nrow = nrow(spectra), ncol = ncol(spectra))
  for (i in 1:nrow(spectra)){
    norm_spectra[i,] <- (spectra[i,] - min(spectra[i,]))/(max(spectra[i,930:1015]) - min(spectra[i,]))
  }
  return(as.data.frame(norm_spectra))
}

