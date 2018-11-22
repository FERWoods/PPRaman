#' Normalising Spectra using SNV
#'
#' Normalises spectra using Standard Normal Variate method
#'
#' @param spectra Spectral like data
#' @returns Normalised spectra using SNV method
#' @export
#'

norm_snv<- function(spectra){

  norm_spec <- matrix(nrow = nrow(spectra), ncol = ncol(spectra))
  for (i in 1:nrow(spectra)){
    norm_spec[i,] <- (spectra[i,] - mean(spectra[i,]))/sd(spectra[i,])
  }

  return(as.data.frame(norm_spec))
}
