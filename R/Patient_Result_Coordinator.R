#' Patient diagnosis
#'
#' Takes individual spectra classification and converts to patient including sens/spec/npv/ppv
#' @param patient_ids Vector containing patient IDs
#' @param clinical_result Vector containing binary clinical result
#' @param raman_result Vector containing binary raman result from classifier
#' @return patient_results, output
#' @export


patient_results <- function(patient_ids, clinical_result, raman_result){
  results <- data.frame("ID" = patient_ids, "Diagnosis" = clinical_result,
                        "Raman" = raman_result)
  results$Diagnosis <- as.numeric(as.character(results$Diagnosis))
  results$Raman <- as.numeric(as.character(results$Raman))
  results$hitormiss <- ifelse(results[,2] == results[,3], 1, 0)
  #Nested if to determine classification outcome
  results$hitmisstype <- ifelse((results$hitormiss == 1 & results[,2] > 0), "TP",
                                ifelse((results$hitormiss == 1 & results[,2] == 0), "TN",
                                       ifelse((results$hitormiss == 0 & results[,2] > 0), "FN",
                                              ifelse((results$hitormiss == 0 & results[,2] == 0), "FP", "HELP"))))
  colnames(results) <- c("ID", "Diagnosis", "Raman", "hitormiss", "hitmisstype")
  result_df <- aggregate(hitmisstype ~ ID, results, paste, collapse = ",") # Stitches all patient reps together
  # Determine patient diagnosis by seeing if the TP or TN wins over FP/FN (i.e. occurs more often for that patient)
  result_df$Patient_Diagnosis <- ifelse(str_count(result_df[,2],"TP") > str_count(result_df[,2], "FN"), "TP",
                                        ifelse(str_count(result_df[,2],"TN") > str_count(result_df[,2], "FP"), "TN",
                                               ifelse(str_count(result_df[,2],"TP") < str_count(result_df[,2], "FN"), "FN",
                                                      ifelse(str_count(result_df[,2],"TN") < str_count(result_df[,2], "FP"), "FP", "NA"))))
  patient_results <- as.data.frame(table(result_df$Patient_Diagnosis))
  #patient_results[2,] <- as.numeric(as.character(patient_results[2,]))
  sensitivity <- patient_results[patient_results$Var1 == "TP", 2]/(patient_results[patient_results$Var1 == "TP",
                                                                                   2] + patient_results[patient_results$Var1 == "FN", 2])
  specificity <-  patient_results[patient_results$Var1 == "TN", 2]/(patient_results[patient_results$Var1 == "TN",
                                                                                    2] + patient_results[patient_results$Var1 == "FP", 2])

  PPV <- patient_results[patient_results$Var1 == "TP", 2]/(patient_results[patient_results$Var1 == "TP",
                                                                           2] + patient_results[patient_results$Var1 == "FP", 2])
  NPV <- patient_results[patient_results$Var1 == "TN", 2]/(patient_results[patient_results$Var1 == "TN",
                                                                           2] + patient_results[patient_results$Var1 == "FN", 2])
  output <- c(sensitivity, specificity, PPV, NPV)
  return(list(result_df, output))
}
