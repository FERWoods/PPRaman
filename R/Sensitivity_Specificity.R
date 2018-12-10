#' Sensitivity Specificity Function
#'
#' To compute sens and spec on a confusion matrix output form RF
#' @param conf_mat Confusion matrix output from random_forest_testing()
#' @return sensitivity and specificity of the model
#' @export

sen_spec <- function(conf_mat){
  sensitivity <- conf_mat[2, 2] / (conf_mat[2, 2] + conf_mat[2, 1]) #TP/TP+FN
  specificity <- conf_mat[1, 1] / (conf_mat[1, 1] + conf_mat[1, 2])
  return(list(sensitivity, specificity))
}

