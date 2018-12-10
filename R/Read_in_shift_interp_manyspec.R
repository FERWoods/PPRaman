#' Reading In, Shifting and Interpolating for many spectra (e.g >100)
#'
#' Reads in spectra, shifts to 1004 and interpolates on the x axis
#' Note currently phenyl peak to 2nd order poly rather than 4th due to error
#'
#'
#'@returns Interpolated spectra ready for pre-processing
#'@export


read_in_shift_interp_manyspec <- function(){
  #user choose folder
  folder <- choose.dir()
  files <- list.files(folder, full.names = TRUE)

  # This reads all the files in using read.table
  inputfiles <- lapply(files, read.table)

  # Extracting date and ID from file name
  files_base <- basename(files)
  flag <- ifelse(substr(files_base, 1, 2) == "Co", 0, 1)

  inx <- regexpr("ID=", files_base) # searching for pattern "ID="
  ID <- substr(files_base, inx +(3), inx +(3+4) ) # 3 to take into account the "ID=", +5 for ID

  # Counts number of repeats for each patient
  rep_count <- as.data.frame(table(ID))

  # Creates a sequence of 1 - nrepeats for labelling
  repeats <- list()
  for (i in 1:nrow(rep_count)){
    repeats[[i]] <- seq(from = 1, to = rep_count[i,2])
  }

  Labels <- cbind(ID, "rep" = unlist(repeats))
  Labels <- cbind(Labels, flag)

  # Means that we have labels eg patient 066.1 066.2 etc
  pat_rep <- paste(Labels[,1], Labels[,2], sep = ".")
  Labels <- cbind(Labels, pat_rep)

  inputfiles <- flip_data(inputfiles)

  #Storing the original wavenumber for creating hyperspec object
  old_wn <- lapply(inputfiles, function(x){ x["V2"] <- NULL; x })

  wn_all <- t(do.call("cbind", old_wn))
  # remove the w/n column from each spectra
  rmv_wn <- lapply(inputfiles, function(x){ x["V1"] <- NULL; x })

  raw_spec <- t(do.call("cbind", rmv_wn))

  #Step sizes for interpolating
  wavenumber <- as.data.frame(seq(611.6, 1717, by = 1.09))

  signal_max <- matrix(nrow = nrow(raw_spec), ncol = 3)
  for (i in 1:nrow(raw_spec)){
    max_i <- which.max(raw_spec[i,300:380]) # taking 1015-300 and 1015-380 due to flipping
    shift_i <- max_i + 299
    signal_max[i,] <- c(shift_i, wn_all[i, shift_i], raw_spec[i, shift_i])
  }

  SP_index <- matrix(nrow = nrow(raw_spec), ncol = 5 ) # using 5 reference points for peak fit
  SP_int <- matrix(nrow = nrow(raw_spec), ncol = 5) # for finding intensity at SP_index
  SP_wn <- matrix(nrow = nrow(raw_spec), ncol = 5) # for wavenumber around the peak
  SP_fit <- list()
  SP_fit_val <- list()

  Shifted_wn <- matrix(nrow = nrow(raw_spec), ncol = ncol(raw_spec))
  Shifted_int <- matrix(nrow = nrow(raw_spec), ncol = ncol(raw_spec))

  ref_peak <- 1004
  for (i in 1:nrow(raw_spec)){
    SP_index[i,] <- c(signal_max[i,1]-2, signal_max[i,1] -1, signal_max[i,1],
                      signal_max[i,1] + 1, signal_max[i,1] + 2)
    SP_int[i,] <- c(raw_spec[i,SP_index[i,1]], raw_spec[i,SP_index[i,2]],
                    raw_spec[i,SP_index[i,3]], raw_spec[i,SP_index[i,4]],
                    raw_spec[i,SP_index[i,5]]) # grabs intensity values around the peak

    SP_wn[i,] <- c(wn_all[i,SP_index[i,1]], wn_all[i,SP_index[i,2]],
                   wn_all[i,SP_index[i,3]], wn_all[i,SP_index[i,4]],
                   wn_all[i,SP_index[i,5]]) # same again for the x axis

    SP_fit[[i]] <- polyfit(SP_wn[i,], SP_int[i,], 2) # fits 2nd order poly to phen peak,
    # 4th order poly causes error with qr.solve -- same problem with claudia, so using 2nd order

    x_ref <- seq(from = SP_wn[i,1], to = SP_wn[i,5], by = 0.1)

    SP_fit_val[[i]] <- polyval(SP_fit[[i]], x_ref)

    max_new_i <- which.max(SP_fit_val[[i]])

    band_pos <- SP_wn[i,1] + 0.1*(max_new_i-1)

    shift <- ref_peak - band_pos

    x_shift <- wn_all[i,] -shift

    Shifted_int[i, ] <- pchip(as.numeric(wn_all[i,]), as.numeric(raw_spec[i,]),
                              as.numeric(x_shift))

    Shifted_wn[i, ] <- x_shift

    interpolated <- matrix(nrow = nrow(raw_spec), ncol = nrow(wavenumber))

  }

  for (i in 1:nrow(raw_spec)){
    interpolated[i,] <- interp1(Shifted_wn[i,], Shifted_int[i,], wavenumber[,1], "linear")

  }

  return(list(interpolated, Labels))
}

