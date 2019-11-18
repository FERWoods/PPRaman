#' Normalise at Phenylalanine peak
#'
#' Normalising function using min/max method at phenylalanine peak
#' @param spectra Spectra in matrix format for normalising
#' @param bin_size Number of wavenumbers in a bin to shift region to normalise to accordingly
#' @return Normalised spectra
#' @export
norm_p <- function(spectra, bin_size = 1){
  norm_spectra <- matrix(nrow = nrow(spectra), ncol = ncol(spectra))
  for (i in 1:nrow(spectra)){
    norm_spectra[i,] <- (spectra[i,] - min(spectra[i,]))/(max(spectra[i,round(300/bin_size):round(400/bin_size)]) - min(spectra[i,]))
  }
  return(as.data.frame(norm_spectra))
}


