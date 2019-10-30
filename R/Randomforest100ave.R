#' Random Forest - Averaging 100
#'
#' Trains model using RF, also includes some QC removing data outside +/- 5 SD
#' @param training_setin This is your training data, with the last column (V1016) to
#' @param testing_setin Testing set with V1016 containing classification independent test set
#' @param cv_setin  Cross validation set
#' @param cv_set_ids Vector with cross val ids
#' @param test_set_ids Vector containing testing set ids
#' @return prob_ave_cv, output_cv, prob_ave_test, output_test
#' @export
#' @import randomForest
#' @import ROCR
random_forest_ave <- function(training_setin, testing_setin, cv_setin,  test_set_ids, cv_set_ids){
  full_train <- rbind(training_setin, cv_setin)
  ## RF generation and averaging over 100
  # CV refers to cross val, whereas the full set includes all spectra to test
  # the GP set
  rf_cv <- list()
  rf_all <- list()
  pred_cv <- list()
  pred_all <- list()
  for (i in 1:100){
    rf_cv[[i]] <- randomForest(V1016 ~.,
                               data = training_setin, importance = TRUE)
    rf_all[[i]] <- randomForest(V1016 ~.,
                                data = full_train, importance = TRUE)
    pred_cv[[i]] <- predict(rf_cv[[i]], cv_setin, type = "prob")
    pred_all[[i]] <- predict(rf_all[[i]], testing_setin, type = "prob")

  }
  # CV set results
  co_ave_cv <- list()
  co_se_cv <- list()
  ca_ave_cv <- list()
  ca_se_cv <- list()
  tmp1_cv <- list()
  tmp2_cv <- list()
  auc_cv <- list()
  for (j in 1:nrow(cv_setin)){
    for (i in 1:100){
      tmp1_cv[[i]] <- pred_cv[[i]][j,1]
      tmp2_cv[[i]] <- pred_cv[[i]][j,2]
      auc_cv[[i]] <- ROCR::performance(prediction(pred_cv[[i]][,2], cv_setin$V1016), "auc")@y.values[[1]]

    }
    co_ave_cv[[j]] <- mean(unlist(tmp1_cv))
    ca_ave_cv[[j]] <- mean(unlist(tmp2_cv))
    co_se_cv[[j]] <- sd(unlist(tmp1_cv))/sqrt(length(unlist(tmp1_cv))) #standard error of the mean
    ca_se_cv[[j]] <- sd(unlist(tmp2_cv))/sqrt(length(unlist(tmp2_cv))) #standard error of the mean
  }

  prob_ave_cv <- data.frame("ID" = cv_set_ids, "Co_prob" = unlist(co_ave_cv), "Ca_prob" = unlist(ca_ave_cv),
                            "Co_SE" = unlist(co_se_cv), "Ca_SE" = unlist(ca_se_cv))
  prob_ave_cv$raman_diagnosis <- ifelse(prob_ave_cv$Co_prob > 0.5, 0, 1)
  prob_ave_cv$clinical_diagnosis <- cv_setin$V1016
  prob_ave_cv$hitmiss <- ifelse(prob_ave_cv$raman_diagnosis == prob_ave_cv$clinical_diagnosis, 1, 0)
  acc_cv <- sum(prob_ave_cv$hitmiss == 1)/length(prob_ave_cv$hitmiss)
  conf_mat_cv <- table(observed = prob_ave_cv$clinical_diagnosis, predicted = prob_ave_cv$raman_diagnosis)
  auc_all_cv <- c(mean(unlist(auc_cv)), sd(unlist(auc_cv))/sqrt(length(unlist(auc_cv))))
  output_cv <- c(acc_cv, sen_spec(conf_mat_cv)[[1]], sen_spec(conf_mat_cv)[[2]], ppv_npv(conf_mat_cv)[[1]],
                 ppv_npv(conf_mat_cv)[[2]], auc_all_cv)

  # Testing set reults
  co_ave_test <- list()
  co_se_test <- list()
  ca_ave_test <- list()
  ca_se_test <- list()
  tmp1_test <- list()
  tmp2_test <- list()
  auc_test <- list()
  for (j in 1:nrow(testing_setin)){
    for (i in 1:100){
      tmp1_test[[i]] <- pred_all[[i]][j,1]
      tmp2_test[[i]] <- pred_all[[i]][j,2]
      auc_test[[i]] <- ROCR::performance(prediction(pred_all[[i]][,2], testing_setin$V1016), "auc")@y.values[[1]]

    }
    co_ave_test[[j]] <- mean(unlist(tmp1_test))
    ca_ave_test[[j]] <- mean(unlist(tmp2_test))
    co_se_test[[j]] <- sd(unlist(tmp1_test))/sqrt(length(unlist(tmp1_test))) #standard error of the mean
    ca_se_test[[j]] <- sd(unlist(tmp2_test))/sqrt(length(unlist(tmp2_test))) #standard error of the mean
  }

  prob_ave_test <- data.frame("ID" = test_set_ids, "Co_prob" = unlist(co_ave_test), "Ca_prob" = unlist(ca_ave_test),
                              "Co_SE" = unlist(co_se_test), "Ca_SE" = unlist(ca_se_test))
  prob_ave_test$raman_diagnosis <- ifelse(prob_ave_test$Co_prob > 0.5, 0, 1)
  prob_ave_test$clinical_diagnosis <- testing_setin$V1016
  prob_ave_test$hitmiss <- ifelse(prob_ave_test$raman_diagnosis == prob_ave_test$clinical_diagnosis, 1, 0)
  acc_test <- sum(prob_ave_test$hitmiss == 1)/length(prob_ave_test$hitmiss)
  conf_mat_test <- table(observed = prob_ave_test$clinical_diagnosis, predicted = prob_ave_test$raman_diagnosis)
  auc_all_test <- c(mean(unlist(auc_test)), sd(unlist(auc_test))/sqrt(length(unlist(auc_test))))

  output_test <- c(acc_test, sen_spec(conf_mat_cv)[[1]], sen_spec(conf_mat_test)[[2]], ppv_npv(conf_mat_test)[[1]],
                   ppv_npv(conf_mat_test)[[2]], auc_all_test)
  return(list(prob_ave_cv, output_cv, prob_ave_test, output_test))
}
