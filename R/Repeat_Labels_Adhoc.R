#' Labels for repeats
#'
#' For use in adhoc type analysis
#'
#' @param label a label for the data, eg 1,1,1,2,2,2
#' @return Returns a data frame contating labels from repeats e.g 1.1, 1.2


reps <- function(label){
  # Counts number of repeats for each patient
  rep_count <- as.data.frame(table(as.data.frame(label)))

  # Creates a sequence of 1 - nrepeats for labelling
  repeats <- list()
  for (i in 1:nrow(rep_count)){
    repeats[[i]] <- seq(from = 1, to = rep_count[i,2])
  }

  Label <- cbind(label, "rep" = unlist(repeats))

  # Means that we have labels eg patient 066.1 066.2 etc
  pat_rep <- paste(Label[,1], Label[,2], sep = ".")
  Labels <- as.data.frame(cbind(Label, pat_rep))

  return(Labels)
}

