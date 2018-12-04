#' Savitsky Golay filter
#'
#' Smoothing using a savitsky golay, parameter is polynomial order
#'
#' @param spectra Spectra data
#' @param poly_order Order of polynomial for filtering
#' @return Filtered spectra
#' @import signal
#' @export

SavGol <- function(spectra, poly_order){
  output <- matrix(ncol = NCOL(spectra), nrow = nrow(spectra))
  for (i in 1:nrow(spectra)){
  output[i,] <- sgolayfilt(spectra[i,], p=poly_order, m=0)
  }
  return(output)
}


