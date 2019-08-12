#' Accuracy calculation for a 2 class model
#'
#' @param predictions Spectra data that has been interpolated already
#' @param ref.labels Reference labels for the spectra
#' @export

calculate.accuracy <- function(predictions, ref.labels) {
  return(length(which(predictions == ref.labels)) / length(ref.labels))
}
