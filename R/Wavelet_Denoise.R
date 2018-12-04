#' Wavelet Denoise
#' Smoothing using the wavelet denoise method
#'
#' @param spectra Spectral object requiring denoising
#' @param filter_len Length of filter to be used with wavelet denoise. must be even and max of N=20
#' @return Denoiseed spectra
#' @references Adapted from PRFFECT package
#' @export


wavelet_denoise <- function(spectra, filter_len){
  Nzeroes<-2^ceiling(log2(NCOL(spectra)))-NCOL(spectra) # for wavelet transform
  output <- matrix(ncol = NCOL(spectra), nrow = nrow(spectra))
  for (i in 1:nrow(spectra)){
    capture.output(suppressMessages(tmp <- as.vector(denoise.dwt(c(spectra[i,], rep(0,times=Nzeroes)),
                                                                 daubcqf(filter_len)$h.0)$xd[1:NCOL(spectra), 1])))
    output[i,] <- tmp
  }
  return(output)
}
