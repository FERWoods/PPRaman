#' Polynomial Baseline Removal
#'
#' Removes baseline contribution using hyperSpec package, more parameters to be added
#'
#'
#' @param spectra Interpolated spectra
#' @param poly.order Polynomial order
#' @return spectra with polynomial baseline removed
#' @export
#'

poly_baseline_removal <- function(spectra, poly.order){
  # Uses hyperSpec polynomial baseline removal
  baseline <- spc.fit.poly.below(spectra, poly.order = poly.order)

  # Remove fitted baseline
  spec_corrected <- spectra - baseline

  return(spec_corrected)

}
