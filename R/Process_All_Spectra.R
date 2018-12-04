#' Processing spectra en masse
#'
#' Runs all spectra in a file through the preprocessing sequence, outputting
#' a data frame of processed spectra.
#'
#' @return A list, with first element the corrected and normalised
#' @return spectra, the other element contrains date + ID information
#' @export

process_all_spec <- function(){

  # Read in all spectra -- list all in one file
  spectral_info <- read_in_shift_interp()
  spectra <- spectral_info[[1]]
  supplementary <- spectral_info[[2]]

  # Run RCF on data
  bl_rmv <- RCF(spectra, nrow(spectra) , "pchip", 150)

  # Normalise to phenylalanine peak
  norm_spec <- norm_p(t(bl_rmv))

  return(list(norm_spec, supplementary))
}
