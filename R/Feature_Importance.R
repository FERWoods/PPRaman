#' Feature Selection for RF
#'
#' Short function to enable rapid extraction of random forest importnace for features
#' @param rf_model Input a random forest model
#' @return rf_imp
#' @export
#' @import randomForest
#' @import ROCR


feature_selection <- function(rf_model){
  rf_imp <- as.data.frame(importance(rf_model, type=2)) # type 2 = MeanGiniDecrease
  rf_imp <- cbind(Feature = rownames(rf_imp), rf_imp)
  rf_imp <- rf_imp[order(rf_imp$MeanDecreaseGini),] # order
  return(rf_imp)
}
