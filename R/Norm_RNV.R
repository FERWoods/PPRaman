#' Normalising Spectra using RNV
#'
#' Normalises spectra using Robust Normal Variate method - RNV was introduced to solve
#' a closure problem - in statistics this is
#'
#' @param spectra Spectral like data
#' @return Normalised spectra using RNV method
#' @export
#'

norm_rnv<- function(spectra){

  norm_spec <- matrix(nrow = nrow(spectra), ncol = ncol(spectra))
  for (i in 1:nrow(spectra)){
    norm_spec[i,] <- (spectra[i,] - mean(spectra[i,]))/sd(spectra[i,])
  }

  return(as.data.frame(norm_spec))
}
