#' Principle Component Analysis
#'
#' Computes PCA on spectra data
#' @param spectra Processed spectral data
#' @return PC loadings and scores in format PC#_score (e.g. pc1_score)
#' @return the full PCA object in the third list object
#' @export

PCA <- function(spectra){

  #Uses base r to complete PCA
  PC <- prcomp(spectra)

  # Assign new variables for each pc loading and put in list
  loadings <- list()
  for (i in 1:ncol(PC$rotation)){
    assign(paste0("pc", i, "_loading"), PC$rotation[,i])
  }

  # Assign new variables for each pc score and put in list
  scores <- list()
  for (i in 1:ncol(PC$x)){
    scores[[i]] <- assign(paste0("pc", i, "_score"), PC$x[,i])
  }

  return(list(scores, loadings, PC))
}
