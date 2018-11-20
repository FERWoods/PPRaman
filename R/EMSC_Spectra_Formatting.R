#' EMSC Spectral Formatting
#'
#' Formats spectra from hyperspec to something consistent with EMSC pkg
#' @param raw_spec Spectra data that has been interpolated already
#' @param n_pat Number of patients - user input
#' @param n_reps Number of repeats for each patient
#' @return Formatted spectra for EMSC
#' @export

format_spectra <- function(raw_spec, n_pat, n_reps){

  #Pre-processing -- need to hammer data into shape to fit the package requirements
  Spectra <- as.data.frame(raw_spec)

  # Gives the steps used in interpolation -- start at 660.65 if
  colnames(Spectra) <- seq(611.6, 1717, by = 1.09)
  Spectra <- cbind(seq(1, 6, by = 1), Spectra)

  # For classifying
  replicates <- unlist(lapply( 1:inputs[[1]], rep, each = inputs[[2]]))
  Spectra <- cbind(Spectra, replicates)

  # Converts type to "AsIs" using I() to work with EMSC
  Spectra_AsIs <- cbind(Spectra[1], Raman = I(as.matrix(Spectra[-1])))
  Spectra_AsIs <- cbind(Spectra_AsIs, replicates = Spectra$replicates)

  #Return formatted spectra
  return(Spectra_AsIs)
}
