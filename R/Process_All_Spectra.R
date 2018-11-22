#' Processing spectra en masse
#'
#' Runs all spectra in a file through the preprocessing sequence, outputting
#' a data frame of processed spectra.
#'
#' @return Processed spectra
#' @export

process_all_spec <- function(){
  # Read in all spectra -- list all in one file
  spectra <- read_interp_spectra()

  # Run RCF on data
  bl_rmv <- RCF(spectra@data$spc, nrow(spectra@data$spc) , "pchip", 150)

  # Normalise to phenylalanine peak
  norm_spec <- norm_p(t(bl_rmv))

  return(norm_spec)
}

