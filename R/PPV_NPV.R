#' NPV PPV Function
#'
#' To compute sens and spec on a confusion matrix output form RF
#' @param conf_mat Confusion matrix output from random_forest_testing()
#' @return PPV NPV of the model
#' @export

ppv_npv <- function(conf_mat){
  ppv <- conf_mat[2, 2] / (conf_mat[2, 2] + conf_mat[1, 2]) #TP/TP+FP
  npv <- conf_mat[1, 1] / (conf_mat[1, 1] + conf_mat[2, 1]) #TN/TN+FN
  return(list(ppv, npv))
}

