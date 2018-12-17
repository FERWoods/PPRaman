#' Random Forest Testing
#'
#' For use with random_forest_training to test models
#'
#' @param rf_model Complete RF model
#' @param rf_model_600 RF model using the top 600 features
#' @param test_set Data to test, last column should be the 1/0 cancer/control flag
#' @param top600 Top 600 features as selected by RF model
#' @return rf_table, rf_table_600 (confusions matrices), rf_perf, rf_perf_600 (ROC curves for each)
#' @export
#' @import randomForest
#' @import ROCR

random_forest_testing <- function(rf_model, rf_model_600, test_set, top600){
  test_set_600 <- test_set[,top600] # Subset training set for 600 top features
  test_set_600$V1016 <- test_set$V1016
  # Predicts using models
  rf_pred <- predict(rf_model, test_set)
  rf_pred_600 <- predict(rf_model_600, test_set_600)
  rf_prob <- predict(rf_model, test_set,type = "prob") # predict (output class probability)
  rf_prob_600 <- predict(rf_model_600, test_set_600,type = "prob") # predict (output class probability)

  # Confusion Matrices
  rf_table <- table(observed = test_set[, ncol(test_set)], predicted = rf_pred)
  rf_table_600 <- table(observed = test_set_600[,ncol(test_set_600)], predicted = rf_pred_600)

  # ROC curves
  rf_perf<-ROCR::performance(prediction(rf_prob[,2],test_set[, ncol(test_set)]),'tpr','fpr')
  rf_perf_600<-ROCR::performance(prediction(rf_prob_600[,2],test_set_600[,ncol(test_set_600)]),'tpr','fpr')

  # sen_spec
  rf_sen_spec <- ROCR::performance(prediction(rf_prob[,2],test_set[, ncol(test_set)]),'sens','spec')
  rf_sen_spec_600 <- ROCR::performance(prediction(rf_prob_600[,2],test_set_600[,ncol(test_set_600)]),'sens','spec')

  # Area under curve
  rf_auc <- ROCR::performance(prediction(rf_prob[,2],test_set[, ncol(test_set)]),'auc')
  rf_auc_600 <- ROCR::performance(prediction(rf_prob_600[,2],test_set_600[,ncol(test_set_600)]),"auc")

  # PPV NPV
  rf_val <- ROCR::performance(prediction(rf_prob[,2],test_set[, ncol(test_set)]),'ppv', 'npv')
  rf_val_600 <- ROCR::performance(prediction(rf_prob_600[,2],test_set_600[,ncol(test_set_600)]),'ppv', 'npv')
  return(list(rf_table, rf_table_600, rf_perf, rf_perf_600, rf_auc, rf_auc_600, rf_sen_spec, rf_sen_spec_600,
              rf_val, rf_val_600))
}
