#' Standard Error
#'
#' @param x list of values to compute SE
#' @return SE
#' @export

std <- function(x){
  sd(x)/sqrt(length(x))
}
