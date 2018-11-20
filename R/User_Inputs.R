#' User Inputs
#'
#' Gathers user inputs on patients/ sample/ repeat numbers
#' @return Numerics for other functions and plotting
#' @export
#'

user_inputs <- function(){
  #user inputs figures
  n_pat <- readline("Number of Patients in sample? Enter into console...")
  n_rep <- readline("Number of Repeats in sample? Enter into console...")
  samp_num <- readline("Total number of samples? Enter into console...")

  #Converting from str to num
  n_pat <- as.numeric(n_pat)
  n_rep <- as.numeric(n_rep)
  samp_num <- as.numeric(samp_num)

  #return numbers as list
  return(list(n_pat, n_rep, samp_num))
}
