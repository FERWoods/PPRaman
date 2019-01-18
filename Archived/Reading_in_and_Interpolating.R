#' Reading in raw spectra and interpolating - ARCHIVED
#'
#' Reads in mutiple raw spectra and interpolates using hyperSpec pkg
#' @return Interpolated spectra in hyperSpec format
#' @export
#' @import hyperSpec

read_interp_spectra <- function(){
  #user chooses
  files <- choose.files()

  # This reads all the files in using read.table
  inputfiles <- lapply(files, read.table)

  inputfiles <- flip_data(inputfiles)
  # Extracting date and ID from file name
  files_base <- basename(files)
  date <- substr(files_base, 1, 6)
  ID <- substr(files_base, 11, 13)
  supp_info <- cbind(date, ID)

  # Counts number of repeats for each patient
  rep_count <- as.data.frame(table(supp_info[,2]))

  # Creates a sequence of 1 - nrepeats for labelling
  repeats <- list()
  for (i in 1:nrow(rep_count)){
    repeats[[i]] <- seq(from = 1, to = rep_count[i,2])
  }

  Labels <- cbind(supp_info, "rep" = unlist(repeats))

  # Means that we have labels eg patient 066.1 066.2 etc
  pat_rep <- paste(Labels[,2], Labels[,3], sep = ".")
  Labels <- cbind(Labels, pat_rep)

  #Storing the original wavenumber for creating hyperspec object
  old_wn <- lapply(inputfiles, function(x){ x["V2"] <- NULL; x })

  # remove the w/n column from each spectra
  rmv_wn <- lapply(inputfiles, function(x){ x["V1"] <- NULL; x })

  raw_spec <- t(do.call("cbind", rmv_wn))

  #Step sizes for interpolating
  waveno <- seq(611.6, 1717, by = 1.09)


  #Converts raw spectra to hyperspec object to work with the package
  raw_hyperSpec <- as.hyperSpec(raw_spec, wl = waveno)
  chondro@wavelength

  #Add labels

  raw_hyperSpec@label$.wavelength <- expression(Wavenumber (cm^-1))
  raw_hyperSpec@label$spc <- expression(Intensity (a.u))

  #Interpolating the data using loess
  int_spec <- spc.loess(raw_hyperSpec, waveno)

  return(list(int_spec, Labels))
}


