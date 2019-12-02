#' Random Forest - Averaging 100
#'
#' Trains model using RF, also includes some QC removing data outside +/- 5 SD
#' @param training_setin This is your training data, with the last column (V1016) to
#' @param testing_setin Testing set with V1016 containing classification independent test set
#' @param test_set_ids Vector containing testing set ids
#' @param training_set_ids Vector containing training set IDs
#' @return prob_ave_cv, output_cv, prob_ave_test, output_test
#' @export
#' @import randomForest
#' @import ROCR
#' @import purrr
#'

#training_setin <- rbind(processed_train_spec, processed_test_spec)
#training_set_ids <- c(as.character(training$ID), as.character(testing$ID))
#testing_setin <- processed_gp_spec
#test_set_ids <- raw_gp_df$ID
random_forest_ave <- function(training_setin, testing_setin, training_set_ids, test_set_ids){
  ## RF generation and averaging over 100
  # CV refers to cross val, whereas the full set includes all spectra to test
  # the GP set
  # preamble
  training_setin <- cbind(training_setin, training_set_ids)
  training_setin$V1016 <- factor(training_setin$V1016)
  colnames(training_setin)[ncol(training_setin)] <- "ID"
  testing_setin <- cbind(testing_setin, test_set_ids)
  testing_setin$V1016 <- factor(testing_setin$V1016)
  colnames(testing_setin)[ncol(testing_setin)] <- "ID"
  ############################################### LOOCV #############################################


  ID <- unique(training_setin$ID)
  loo_ca <- list()
  for(i in ID){
    rf <- randomForest(V1016~., data=training_setin[training_setin$ID!=i,1:(ncol(training_setin)-1)])
    loo_ca[[i]] <- predict(rf, newdata=training_setin[training_setin$ID==i,1:(ncol(training_setin)-1)],
                        type = "prob")[,2]
  }
  res <- data.frame("ID" = training_set_ids, "Raman_Ca_Prob" = unlist(loo_ca),
                    "Clinical_Diagnosis" = as.numeric(as.character(training_setin$V1016)),
                    "Raman_Diagnosis" = ifelse(unlist(loo_ca) >= 0.5, 1, 0))

  # Spectra by spectra results
  loocv_conf_mat_spec <- table(observed = res$Clinical_Diagnosis, predicted = res$Raman_Diagnosis)
  loocv_senspec_spec <- sen_spec(loocv_conf_mat_spec)
  loocv_ppvnpv_spec <- ppv_npv(loocv_conf_mat_spec)
  loocv_acc_spec <- calculate.accuracy(res$Raman_Diagnosis, res$Clinical_Diagnosis)

  #Aggregating on patient IDs and averaging resultant probability
  loocv_res <- as.data.frame(aggregate(Raman_Ca_Prob~ID, res, mean))
  loocv_res$Raman_Diagnosis <- ifelse(loocv_res$Raman_Ca_Prob >= 0.5, 1, 0)
  loocv_res$SD <- aggregate(Raman_Ca_Prob~ID, res, sd)[,2]
  loocv_res$Clinical_Diagnosis <- aggregate(Clinical_Diagnosis~ID, res, mean)[,2]
  loocv_conf_mat_pat <- table(observed = loocv_res$Clinical_Diagnosis, predicted = loocv_res$Raman_Diagnosis)
  loocv_senspec_pat <- sen_spec(loocv_conf_mat_pat)
  loocv_ppvnpv_pat <- ppv_npv(loocv_conf_mat_pat)
  loocv_acc_pat <- calculate.accuracy(loocv_res$Raman_Diagnosis, loocv_res$Clinical_Diagnosis)

  # Training_set variance
  mean_mat <- matrix(nrow=length(ID), ncol = (ncol(training_setin)-2))
  sd_mat <- matrix(nrow=length(ID), ncol = (ncol(training_setin)-2))
  for (i in 1:length(ID)){
    mean_mat[i,] <- apply(training_setin[training_setin$ID == ID[i],1:(ncol(training_setin)-2)], 2, mean) # minus 2 to rmv id and binary
    sd_mat[i,] <- apply(training_setin[training_setin$ID == ID[i],1:(ncol(training_setin)-2)], 2, sd)
  }

  train_spec_sd <- rowMeans(sd_mat) # SCOPE TO GO BACK AND RETEST THESE TO PULL ALL THE WN FOR TOP 10 SAY
  loocv_res$specSD <- train_spec_sd

  loocv_out <- cbind(loocv_senspec_spec[[1]], loocv_senspec_spec[[2]], loocv_ppvnpv_spec[[1]], loocv_ppvnpv_spec[[2]],
                     loocv_acc_spec[1], loocv_senspec_pat[[1]], loocv_senspec_pat[[2]], loocv_ppvnpv_pat[[1]],
                     loocv_ppvnpv_pat[[2]], loocv_acc_pat[1], mean(train_spec_sd), mean(loocv_res$SD))

  ############################################### LOOCV END ##############################################

  ############################################### Testing Set ##############################################

  rf_all <- list()
  ca_all <- list()
  rf_imp_tmp <- list()
  rf_imp <- list()
  for (i in 1:100){
    rf_all[[i]] <- randomForest(V1016 ~.,
                                data = training_setin[,1:(ncol(training_setin)-1)], importance = TRUE)
    ca_all[[i]] <- predict(rf_all[[i]], testing_setin[,1:(ncol(testing_setin)-2)], type = "prob")[,2]
    # Find top 600 features
    rf_imp_tmp[[i]] <- as.data.frame(importance(rf_all[[i]], type=2)) # type 2 = MeanGiniDecrease
    rf_imp[[i]] <- cbind(Feature = rownames(rf_imp_tmp[[i]]), rf_imp_tmp[[i]])

  }
  imp_all <- lapply(rf_imp, setNames, c("Feature", "MeanDecreaseGini"))
  imp_all_comb <- reduce(imp_all, full_join, by = "Feature")
  imp_res <- cbind(as.character(imp_all_comb[,1]),rowMeans(imp_all_comb[,2:101]), apply(imp_all_comb[,2:101], 1, sd))

  # Testing set reults
  ca_ave_test <- list()
  ca_se_test <- list()
  tmp1_test <- list()
  auc_test <- list()
  for (j in 1:nrow(testing_setin)){
    for (i in 1:100){
      tmp1_test[[i]] <- ca_all[[i]][j]

      auc_test[[i]] <- ROCR::performance(prediction(ca_all[[i]], testing_setin$V1016), "auc")@y.values[[1]]

    }
    ca_ave_test[[j]] <- mean(unlist(tmp1_test))
    ca_se_test[[j]] <- sd(unlist(tmp1_test))/sqrt(length(unlist(tmp1_test))) #standard error of the mean
  }

  res_test <- data.frame("ID" = test_set_ids, "Raman_Ca_Prob" = unlist(ca_ave_test),
                         "Raman_Ca_Prob_SD" = unlist(ca_se_test),
                    "Clinical_Diagnosis" = as.numeric(as.character(testing_setin$V1016)),
                    "Raman_Diagnosis" = ifelse(unlist(ca_ave_test) >= 0.5, 1, 0))

  auc_ave_test <- mean(unlist(auc_test))
  auc_se_test <- sd(unlist(auc_test))/sqrt(length(unlist(auc_test))) #standard error of the mean
  # Spectra by spectra results
  test_conf_mat_spec <- table(observed = res_test$Clinical_Diagnosis, predicted = res_test$Raman_Diagnosis)
  test_senspec_spec <- sen_spec(test_conf_mat_spec)
  test_ppvnpv_spec <- ppv_npv(test_conf_mat_spec)
  test_acc_spec <- calculate.accuracy(res_test$Raman_Diagnosis, res_test$Clinical_Diagnosis)

  #Aggregating on patient IDs and averaging resultant probability
  test_res <- as.data.frame(aggregate(Raman_Ca_Prob~ID, res_test, mean))
  test_res$Raman_Diagnosis <- ifelse(test_res$Raman_Ca_Prob >= 0.5, 1, 0)
  test_res$SD <- aggregate(Raman_Ca_Prob~ID, res_test, sd)[,2]
  test_res$Clinical_Diagnosis <- aggregate(Clinical_Diagnosis~ID, res_test, mean)[,2]
  test_conf_mat_pat <- table(observed = test_res$Clinical_Diagnosis, predicted = test_res$Raman_Diagnosis)
  test_senspec_pat <- sen_spec(test_conf_mat_pat)
  test_ppvnpv_pat <- ppv_npv(test_conf_mat_pat)
  test_acc_pat <- calculate.accuracy(test_res$Raman_Diagnosis, test_res$Clinical_Diagnosis)

  # Training_set variance
  test_ID <- unique(test_set_ids)
  test_mean_mat <- matrix(nrow=length(test_ID), ncol = (ncol(testing_setin)-2))
  test_sd_mat <- matrix(nrow=length(test_ID), ncol = (ncol(testing_setin)-2))
  for (i in 1:length(test_ID)){
    test_mean_mat[i,] <- apply(testing_setin[testing_setin$ID == test_ID[i],1:(ncol(testing_setin)-2)], 2, mean)
    test_sd_mat[i,] <- apply(testing_setin[testing_setin$ID == test_ID[i],1:(ncol(testing_setin)-2)], 2, sd)
  }

  test_spec_sd <- rowMeans(test_sd_mat)
  test_res$specSD <- test_spec_sd

  test_out <- cbind(test_senspec_spec[[1]], test_senspec_spec[[2]], test_ppvnpv_spec[[1]], test_ppvnpv_spec[[2]],
                    test_acc_spec[1], test_senspec_pat[[1]], test_senspec_pat[[2]], test_ppvnpv_pat[[1]],
                    test_ppvnpv_pat[[2]], test_acc_pat[1], mean(test_spec_sd), mean(test_res$SD), auc_ave_test, auc_se_test)

  ############################################### Testing Set End ##############################################

  return(list(loocv_out, res, loocv_res, test_out, res_test, test_res, imp_res))
}
