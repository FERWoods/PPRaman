#' PC Plots
#'
#' Loadings and score plots ready to use
#' @param PCA_output Output from using PCA() function
#' @return Data frame of PC1 and PC2 and their labels for plotting
#' @export

pc_plots <- function(PCA_output){
  # Extract PC1 and PC2 scores from PCA output
  PC1_PC2_Scores <- cbind(as.data.frame(unlist(PCA_output[[1]][1])),
                          as.data.frame(unlist(PCA_output[[1]][2])))
  #Combine with patient label for plotting
  labels <- as.factor(unlist(lapply( 1:inputs[[1]], rep, each = inputs[[2]])))
  PC1_PC2_Scores <- cbind(labels, PC1_PC2_Scores)
  colnames(PC1_PC2_Scores) <- c("Patient","PC1", "PC2")

  #Returns dataframe ready for plotting
  return(PC1_PC2_Scores)
}
