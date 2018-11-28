#' Normalising Spectra using Vector Normalising
#'
#'
#' @param spectra Spectral like data
#' @return Normalised spectra using vector normalising
#' @export
#'

norm_vec<- function(spectra){

  norm_spec <- matrix(nrow = nrow(spectra), ncol = ncol(spectra))
  for (i in 1:nrow(spectra)){
    norm_spec[i,] <- spectra[i,]/ norm(spectra[i,], type = "2")
  }

  return(as.data.frame(norm_spec))
}
