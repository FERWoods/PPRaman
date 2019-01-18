#' HPC READ IN AND INTERPOLATE SPECTRA
#'
#' Reads in spectra, shifts to 1004 and interpolates on the x axis
#' Note currently phenyl peak to 2nd order poly rather than 4th due to error
#'
#'@returns Interpolated spectra ready for pre-processing
#'@export


hpc_read_in_interp <- function(spectra){
  # Ensure spectra
  inputfiles <- flip_data(spectra)

  #Storing the original wavenumber for creating hyperspec object
  old_wn <- lapply(inputfiles, function(x){ x["V2"] <- NULL; x })

  wn_all <- t(do.call("cbind", old_wn))
  # remove the w/n column from each spectra
  rmv_wn <- lapply(inputfiles, function(x){ x["V1"] <- NULL; x })

  raw_spec <- t(do.call("cbind", rmv_wn))

  #Step sizes for interpolating
  wavenumber <- as.data.frame(seq(611.6, 1717, by = 1.09))

  signal_max <- matrix(nrow = 6, ncol = 3)
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

    Shifted_int[i, ] <- pchip(wn_all[i,], raw_spec[i,], x_shift)
    Shifted_wn[i, ] <- x_shift

    interpolated <- matrix(nrow = nrow(raw_spec), ncol = nrow(wavenumber))

  }
  for (i in 1:nrow(raw_spec)){
    interpolated[i,] <- interp1(Shifted_wn[i,], Shifted_int[i,], wavenumber[,1], "linear")

  }

  return(list(interpolated, Labels))
}
