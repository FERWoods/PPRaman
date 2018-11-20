#' Normalise at Phenylalanine peak
#'
#' Normalising function using min/max method at phenylalanine peak
#' @param spectra Spectra in matrix format for normalising
#' @return Normalised spectra
#' @export
norm_p <- function(spectra){
  norm_spectra <- matrix(nrow = nrow(spectra), ncol = ncol(spectra))
  for (i in 1:nrow(spectra)){
    norm_spectra[i,] <- (spectra[i,] - min(spectra[i,]))/(max(spectra[i,300:400]) - min(spectra[i,]))
  }
  return(as.data.frame(norm_spectra))
}


