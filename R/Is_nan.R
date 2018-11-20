#' Replace NaN with 0
#' Function for replacing NaN's with zeros in the entire df
#' @param x Data frame / structure
#' @return Data frame with all NaN's replaced with 0's
#' @export

is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan))

}
