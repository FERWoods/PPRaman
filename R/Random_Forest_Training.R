#' Random Forest Training
#'
#' Trains model using RF, also includes some QC removing data outside +/- 5 SD
#' @param training_setin This is your training data, with the last column (V1016) to
#' @param training_setin contain 1 or 0 for cancer/control respectively
#' @param rmv_out Choice to remove outliers > 5SD from data points (and if more than 5 points at >5SDs)
#' @return rf_model, rf_model_600 (uses top 600 features), top600 (top 600 features)
#' @export
#' @import randomForest
#' @import ROCR

random_forest_training <- function(training_setin, rmv_out){
  if (rmv_out == 1){
    # Ensure cancer flag is factor - assumes last column in training set is cancer flag
    training_setin[,ncol(training_setin)]<-as.factor(training_setin[,ncol(training_setin)])

    # Function finds the mean and sd of x and outputs logical if x is outside 5 SDs
    remove_outliers <- function(x, limit = 5) {
      mn <- mean(x, na.rm = TRUE)
      out <- limit * sd(x, na.rm = TRUE)
      x < (mn - out) | x > (mn + out)
    }

    # Applies across each column and outputs true or false
    outliers <- apply(training_setin[,-ncol(training_setin)], 2, FUN = remove_outliers)

    # Finds index of any row with outliers
    outliers2 <- which(outliers == TRUE, arr.ind = TRUE)
    tally <- as.data.frame(table(outliers2[,1]))

    row_out_ind <- subset(tally, tally$Freq >= 5)
    index <- as.numeric(as.character(row_out_ind$Var1))

    # remove rows with outliers
    training_rmv_outliers <- training_setin[-index,]

    # Applying RF
    rf_model <- randomForest(V1016 ~.,
                             data = training_rmv_outliers, importance = TRUE)

    # Find top 600 features
    rf_imp <- as.data.frame(importance(rf_model, type=2)) # type 2 = MeanGiniDecrease
    rf_imp <- cbind(Feature = rownames(rf_imp), rf_imp)
    rf_imp <- rf_imp[order(rf_imp$MeanDecreaseGini),] # order
    top600 <- rf_imp$Feature[1:600]

    new_train_600 <- as.data.frame(training_rmv_outliers[, top600])
    new_train_600$V1016 <- training_rmv_outliers$V1016

    # Retraining using top 600 features
    rf_model_600 <- randomForest(V1016 ~.,data = new_train_600,importance=TRUE) #train
    return(list(rf_model, rf_model_600, top600, row_out_ind, training_rmv_outliers))

  } else if(rmv_out == 0){

    # Ensure cancer flag is factor - assumes last column in training set is cancer flag
    training_setin[,ncol(training_setin)]<-as.factor(training_setin[,ncol(training_setin)])

    # Applying RF
    rf_model <- randomForest(V1016 ~.,
                             data = training_setin, importance = TRUE)

    # Find top 600 features
    rf_imp <- as.data.frame(importance(rf_model, type=2)) # type 2 = MeanGiniDecrease
    rf_imp <- cbind(Feature = rownames(rf_imp), rf_imp)
    rf_imp <- rf_imp[order(rf_imp$MeanDecreaseGini),] # order
    top600 <- rf_imp$Feature[1:600]

    new_train_600 <- as.data.frame(training_setin[, top600])
    new_train_600$V1016 <- training_setin$V1016

    # Retraining using top 600 features
    rf_model_600 <- randomForest(V1016 ~.,data = new_train_600,importance=TRUE) #train
    return(list(rf_model, rf_model_600, top600, row_out_ind, training_setin))

  }
}
