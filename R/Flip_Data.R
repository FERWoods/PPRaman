#' Flip data
#'
#' Flips the data so it runs from low to high
#' @param files Raw input files
#' @return The input files flipped from low to high
#' @export

flip_data <- function(files){
  files_flipped <- list()
  for(i in 1:length(files)){
    if (files[[i]][1,1] > files[[i]][nrow(files[[i]]), 1]){
      files_flipped[[i]] <- files[[i]][nrow(files[[i]]):1,]
    }
    else {
      files_flipped[[i]] <- files[[i]]
    }
  }
  return(files_flipped)
}

