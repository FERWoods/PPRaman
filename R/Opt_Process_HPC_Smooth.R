#' Options for pre-processing spectra - developing version
#'
#' @param spectra Spectra data that has been interpolated already
#' @param norm_meth Choice of normalisation
#' @param rm_bl Choice of background removal
#' @param RCF_rad RCF radius selection
#' @param poly_bl_order Choice of polynomial order for baseline removal
#' @param deriv_order derivative order choice for baseline correction
#' @param filter_length filter length (window) for savitsky golay smoothing
#' @param poly_smooth_order polynomial order for smoothing
#' @param smoothing_opt Binary 1 or 0 if you want smoothing or not! This also automatically removes the dud wavenumbers from the fl + 1/2 at either end
#' @param binning_opt Binary 1 or 0 for binning option
#' @param bin_size A numeric value for the desired # of wavenumbers in each bin
#' @param pp_order 1 = Binning (1 must always come first), 2 = smoothing parameters, 3 = baseline correction, 4 = normalisation
#' @return Preprocessed spectra
#' @export

opt_process_hpc_smooth <- function(spectra, norm_meth, rm_bl, RCF_rad, poly_bl_order,
                                   filter_length, deriv_order = 0, poly_smooth_order, smoothing_opt = 0,
                                   binning_opt = 0, bin_size = 1, pp_order = "1234"){
  if(pp_order == "1234"){
    if(binning_opt == 1){
      spectra_binned <- apply(spectra[, 1:1015], 1, binning, bin.size = bin_size)
      wavenumber <- as.data.frame(binning(wavenumber[,1], bin.size = bin_size))
      spectra <- t(spectra_binned)
    }

    if(smoothing_opt == 1){
      smooth_spectra <- apply(spectra, 1, savgol, fl = as.numeric(filter_length),
                              forder = as.numeric(poly_smooth_order), dorder = as.numeric(deriv_order))
      # Step required since we lose some wavenumbers due to filter length
      smooth_tmp <- t(smooth_spectra)
      smooth_spectra <- smooth_tmp[,-c(1:((filter_length + 1)/2),
                                           ((ncol(smooth_tmp) - (1+(filter_length + 1)/2)):ncol(smooth_tmp)))]

      spectra <- smooth_spectra
      wavenumber <- data.frame(wavenumber[-c(1:((filter_length + 1)/2),
                                             ((ncol(smooth_tmp) - (1+(filter_length + 1)/2)):ncol(smooth_tmp))),])
    }

    # Baseline removal selection
    if(rm_bl == "rcf"){
      spectra <- t(apply(spectra, 1, RCF_GENERALISED, radius = as.numeric(RCF_rad), wavenumber = wavenumber))
    } else if(rm_bl == "der"){
      spectra <- spectra # would have already applied the derivative bl removal in the smoothing step
    } else if(rm_bl == "rub"){
      temp <- spc.rubberband(spectra)
      spectra <- spectra - temp
    } else if(rm_bl == "pol"){
      spectra <- t(poly_baseline_removal(as.hyperSpec(spectra), poly.order = poly_bl_order))
      #smooth_replacement <- zeros(m = (3 + 1)/2, n = nrow(spectra))
      #smooth_replacement2 <- cbind(smooth_replacement, spectra, smooth_replacement)
      #spectra <- smooth_replacement2
    } else if(rm_bl == "none"){
      spectra <- spectra
    }

    # Normalise selection
    if(norm_meth == "norm_p"){
      spectra <- norm_p(spectra, bin_size = bin_size)
    } else if(norm_meth == "norm_minmax"){
      spectra <- norm_minmax(spectra)
    } else if(norm_meth == "norm_snv"){
      spectra <- norm_snv(spectra)
    } else if(norm_meth == "norm_vec"){
      spectra <- norm_vec(spectra)
    } else if(norm_meth == "norm_a"){
      spectra <- norm_a(spectra, bin_size = bin_size)
    } else if(norm_meth == "none"){
      spectra <- spectra
    }

    return(as.data.frame(spectra))
  } else if(pp_order == "1324"){
    # Binning
    if(binning_opt == 1){
      spectra_binned <- apply(spectra[, 1:1015], 1, binning, bin.size = bin_size)
      wavenumber <- as.data.frame(binning(wavenumber[,1], bin.size = bin_size))
      spectra <- t(spectra_binned)
    }
    # Baseline removal selection
    if(rm_bl == "rcf"){
      spectra <- t(apply(spectra, 1, RCF_GENERALISED, radius = as.numeric(RCF_rad), wavenumber = wavenumber))
    } else if(rm_bl == "der"){
      spectra <- spectra # would have already applied the derivative bl removal in the smoothing step
    } else if(rm_bl == "rub"){
      temp <- spc.rubberband(spectra)
      spectra <- spectra - temp
    } else if(rm_bl == "pol"){
      spectra <- t(poly_baseline_removal(as.hyperSpec(spectra), poly.order = poly_bl_order))
      #smooth_replacement <- zeros(m = (3 + 1)/2, n = nrow(spectra))
      #smooth_replacement2 <- cbind(smooth_replacement, spectra, smooth_replacement)
      #spectra <- smooth_replacement2
    } else if(rm_bl == "none"){
      spectra <- spectra
    }

    if(smoothing_opt == 1){
      smooth_spectra <- apply(spectra, 1, savgol, fl = as.numeric(filter_length),
                              forder = as.numeric(poly_smooth_order), dorder = as.numeric(deriv_order))
      # Step required since we lose some wavenumbers due to filter length
      smooth_tmp <- t(smooth_spectra)
      smooth_spectra <- smooth_tmp[,-c(1:((filter_length + 1)/2),
                                       ((ncol(smooth_tmp) - (1+(filter_length + 1)/2)):ncol(smooth_tmp)))]

      spectra <- smooth_spectra
      wavenumber <- data.frame(wavenumber[-c(1:((filter_length + 1)/2),
                                             ((ncol(smooth_tmp) - (1+(filter_length + 1)/2)):ncol(smooth_tmp))),])
    }


    # Normalise selection
    if(norm_meth == "norm_p"){
      spectra <- norm_p(spectra, bin_size = bin_size)
    } else if(norm_meth == "norm_minmax"){
      spectra <- norm_minmax(spectra)
    } else if(norm_meth == "norm_snv"){
      spectra <- norm_snv(spectra)
    } else if(norm_meth == "norm_vec"){
      spectra <- norm_vec(spectra)
    } else if(norm_meth == "norm_a"){
      spectra <- norm_a(spectra, bin_size = bin_size)
    } else if(norm_meth == "none"){
      spectra <- spectra
    }

    return(as.data.frame(spectra))

  } else if(pp_order == "1243"){
    #Binning
    if(binning_opt == 1){
      spectra_binned <- apply(spectra[, 1:1015], 1, binning, bin.size = bin_size)
      wavenumber <- as.data.frame(binning(wavenumber[,1], bin.size = bin_size))
      spectra <- t(spectra_binned)
    }

    if(smoothing_opt == 1){
      smooth_spectra <- apply(spectra, 1, savgol, fl = as.numeric(filter_length),
                              forder = as.numeric(poly_smooth_order), dorder = as.numeric(deriv_order))
      # Step required since we lose some wavenumbers due to filter length
      smooth_tmp <- t(smooth_spectra)
      smooth_spectra <- smooth_tmp[,-c(1:((filter_length + 1)/2),
                                       ((ncol(smooth_tmp) - (1+(filter_length + 1)/2)):ncol(smooth_tmp)))]

      spectra <- smooth_spectra
      wavenumber <- data.frame(wavenumber[-c(1:((filter_length + 1)/2),
                                             ((ncol(smooth_tmp) - (1+(filter_length + 1)/2)):ncol(smooth_tmp))),])
    }

    # Normalise selection
    if(norm_meth == "norm_p"){
      spectra <- norm_p(spectra, bin_size = bin_size)
    } else if(norm_meth == "norm_minmax"){
      spectra <- norm_minmax(spectra)
    } else if(norm_meth == "norm_snv"){
      spectra <- norm_snv(spectra)
    } else if(norm_meth == "norm_vec"){
      spectra <- norm_vec(spectra)
    } else if(norm_meth == "norm_a"){
      spectra <- norm_a(spectra, bin_size = bin_size)
    } else if(norm_meth == "none"){
      spectra <- spectra
    }

    # Baseline removal selection
    if(rm_bl == "rcf"){
      spectra <- t(apply(spectra, 1, RCF_GENERALISED, radius = as.numeric(RCF_rad), wavenumber = wavenumber))
    } else if(rm_bl == "der"){
      spectra <- spectra # would have already applied the derivative bl removal in the smoothing step
    } else if(rm_bl == "rub"){
      temp <- spc.rubberband(spectra)
      spectra <- spectra - temp
    } else if(rm_bl == "pol"){
      spectra <- t(poly_baseline_removal(as.hyperSpec(spectra), poly.order = poly_bl_order))
      #smooth_replacement <- zeros(m = (3 + 1)/2, n = nrow(spectra))
      #smooth_replacement2 <- cbind(smooth_replacement, spectra, smooth_replacement)
      #spectra <- smooth_replacement2
    } else if(rm_bl == "none"){
      spectra <- spectra
    }

    return(as.data.frame(spectra))
  } else if(pp_order == "1342"){
    if(binning_opt == 1){
      spectra_binned <- apply(spectra[, 1:1015], 1, binning, bin.size = bin_size)
      wavenumber <- as.data.frame(binning(wavenumber[,1], bin.size = bin_size))
      spectra <- t(spectra_binned)
    }
    # Baseline removal selection
    if(rm_bl == "rcf"){
      spectra <- t(apply(spectra, 1, RCF_GENERALISED, radius = as.numeric(RCF_rad), wavenumber = wavenumber))
    } else if(rm_bl == "der"){
      spectra <- spectra # would have already applied the derivative bl removal in the smoothing step
    } else if(rm_bl == "rub"){
      temp <- spc.rubberband(spectra)
      spectra <- spectra - temp
    } else if(rm_bl == "pol"){
      spectra <- t(poly_baseline_removal(as.hyperSpec(spectra), poly.order = poly_bl_order))
      #smooth_replacement <- zeros(m = (3 + 1)/2, n = nrow(spectra))
      #smooth_replacement2 <- cbind(smooth_replacement, spectra, smooth_replacement)
      #spectra <- smooth_replacement2
    } else if(rm_bl == "none"){
      spectra <- spectra
    }

    # Normalise selection
    if(norm_meth == "norm_p"){
      spectra <- norm_p(spectra, bin_size = bin_size)
    } else if(norm_meth == "norm_minmax"){
      spectra <- norm_minmax(spectra)
    } else if(norm_meth == "norm_snv"){
      spectra <- norm_snv(spectra)
    } else if(norm_meth == "norm_vec"){
      spectra <- norm_vec(spectra)
    } else if(norm_meth == "norm_a"){
      spectra <- norm_a(spectra, bin_size = bin_size)
    } else if(norm_meth == "none"){
      spectra <- spectra
    }

    if(smoothing_opt == 1){
      smooth_spectra <- apply(spectra, 1, savgol, fl = as.numeric(filter_length),
                              forder = as.numeric(poly_smooth_order), dorder = as.numeric(deriv_order))
      # Step required since we lose some wavenumbers due to filter length
      smooth_tmp <- t(smooth_spectra)
      smooth_spectra <- smooth_tmp[,-c(1:((filter_length + 1)/2),
                                       ((ncol(smooth_tmp) - (1+(filter_length + 1)/2)):ncol(smooth_tmp)))]

      spectra <- smooth_spectra
      wavenumber <- data.frame(wavenumber[-c(1:((filter_length + 1)/2),
                                             ((ncol(smooth_tmp) - (1+(filter_length + 1)/2)):ncol(smooth_tmp))),])
    }

    return(as.data.frame(spectra))

  } else if(pp_order == "1432"){

    if(binning_opt == 1){
      spectra_binned <- apply(spectra[, 1:1015], 1, binning, bin.size = bin_size)
      wavenumber <- as.data.frame(binning(wavenumber[,1], bin.size = bin_size))
      spectra <- t(spectra_binned)
    }

    # Normalise selection
    if(norm_meth == "norm_p"){
      spectra <- norm_p(spectra, bin_size = bin_size)
    } else if(norm_meth == "norm_minmax"){
      spectra <- norm_minmax(spectra)
    } else if(norm_meth == "norm_snv"){
      spectra <- norm_snv(spectra)
    } else if(norm_meth == "norm_vec"){
      spectra <- norm_vec(spectra)
    } else if(norm_meth == "norm_a"){
      spectra <- norm_a(spectra, bin_size = bin_size)
    } else if(norm_meth == "none"){
      spectra <- spectra
    }

    # Baseline removal selection
    if(rm_bl == "rcf"){
      spectra <- t(apply(spectra, 1, RCF_GENERALISED, radius = as.numeric(RCF_rad), wavenumber = wavenumber))
    } else if(rm_bl == "der"){
      spectra <- spectra # would have already applied the derivative bl removal in the smoothing step
    } else if(rm_bl == "rub"){
      temp <- spc.rubberband(spectra)
      spectra <- spectra - temp
    } else if(rm_bl == "pol"){
      spectra <- t(poly_baseline_removal(as.hyperSpec(spectra), poly.order = poly_bl_order))
      #smooth_replacement <- zeros(m = (3 + 1)/2, n = nrow(spectra))
      #smooth_replacement2 <- cbind(smooth_replacement, spectra, smooth_replacement)
      #spectra <- smooth_replacement2
    } else if(rm_bl == "none"){
      spectra <- spectra
    }

    if(smoothing_opt == 1){
      smooth_spectra <- apply(spectra, 1, savgol, fl = as.numeric(filter_length),
                              forder = as.numeric(poly_smooth_order), dorder = as.numeric(deriv_order))
      # Step required since we lose some wavenumbers due to filter length
      smooth_tmp <- t(smooth_spectra)
      smooth_spectra <- smooth_tmp[,-c(1:((filter_length + 1)/2),
                                       ((ncol(smooth_tmp) - (1+(filter_length + 1)/2)):ncol(smooth_tmp)))]

      spectra <- smooth_spectra
      wavenumber <- data.frame(wavenumber[-c(1:((filter_length + 1)/2),
                                             ((ncol(smooth_tmp) - (1+(filter_length + 1)/2)):ncol(smooth_tmp))),])
    }


    return(as.data.frame(spectra))
  } else if(pp_order == "1423"){
    # Binning
    if(binning_opt == 1){
      spectra_binned <- apply(spectra[, 1:1015], 1, binning, bin.size = bin_size)
      wavenumber <- as.data.frame(binning(wavenumber[,1], bin.size = bin_size))
      spectra <- t(spectra_binned)
    }

    # Normalise selection
    if(norm_meth == "norm_p"){
      spectra <- norm_p(spectra, bin_size = bin_size)
    } else if(norm_meth == "norm_minmax"){
      spectra <- norm_minmax(spectra)
    } else if(norm_meth == "norm_snv"){
      spectra <- norm_snv(spectra)
    } else if(norm_meth == "norm_vec"){
      spectra <- norm_vec(spectra)
    } else if(norm_meth == "norm_a"){
      spectra <- norm_a(spectra, bin_size = bin_size)
    } else if(norm_meth == "none"){
      spectra <- spectra
    }
    if(smoothing_opt == 1){
      smooth_spectra <- apply(spectra, 1, savgol, fl = as.numeric(filter_length),
                              forder = as.numeric(poly_smooth_order), dorder = as.numeric(deriv_order))
      # Step required since we lose some wavenumbers due to filter length
      smooth_tmp <- t(smooth_spectra)
      smooth_spectra <- smooth_tmp[,-c(1:((filter_length + 1)/2),
                                       ((ncol(smooth_tmp) - (1+(filter_length + 1)/2)):ncol(smooth_tmp)))]

      spectra <- smooth_spectra
      wavenumber <- data.frame(wavenumber[-c(1:((filter_length + 1)/2),
                                             ((ncol(smooth_tmp) - (1+(filter_length + 1)/2)):ncol(smooth_tmp))),])
    }

    # Baseline removal selection
    if(rm_bl == "rcf"){
      spectra <- t(apply(spectra, 1, RCF_GENERALISED, radius = as.numeric(RCF_rad), wavenumber = wavenumber))
    } else if(rm_bl == "der"){
      spectra <- spectra # would have already applied the derivative bl removal in the smoothing step
    } else if(rm_bl == "rub"){
      temp <- spc.rubberband(spectra)
      spectra <- spectra - temp
    } else if(rm_bl == "pol"){
      spectra <- t(poly_baseline_removal(as.hyperSpec(spectra), poly.order = poly_bl_order))
      #smooth_replacement <- zeros(m = (3 + 1)/2, n = nrow(spectra))
      #smooth_replacement2 <- cbind(smooth_replacement, spectra, smooth_replacement)
      #spectra <- smooth_replacement2
    } else if(rm_bl == "none"){
      spectra <- spectra
    }

    return(as.data.frame(spectra))
  } else if(pp_order == "1324"){
    # Binning
    if(binning_opt == 1){
      spectra_binned <- apply(spectra[, 1:1015], 1, binning, bin.size = bin_size)
      wavenumber <- as.data.frame(binning(wavenumber[,1], bin.size = bin_size))
      spectra <- t(spectra_binned)
    }
    # Baseline removal selection
    if(rm_bl == "rcf"){
      spectra <- t(apply(spectra, 1, RCF_GENERALISED, radius = as.numeric(RCF_rad), wavenumber = wavenumber))
    } else if(rm_bl == "der"){
      spectra <- spectra # would have already applied the derivative bl removal in the smoothing step
    } else if(rm_bl == "rub"){
      temp <- spc.rubberband(spectra)
      spectra <- spectra - temp
    } else if(rm_bl == "pol"){
      spectra <- t(poly_baseline_removal(as.hyperSpec(spectra), poly.order = poly_bl_order))
      #smooth_replacement <- zeros(m = (3 + 1)/2, n = nrow(spectra))
      #smooth_replacement2 <- cbind(smooth_replacement, spectra, smooth_replacement)
      #spectra <- smooth_replacement2
    } else if(rm_bl == "none"){
      spectra <- spectra
    }

    if(smoothing_opt == 1){
      smooth_spectra <- apply(spectra, 1, savgol, fl = as.numeric(filter_length),
                              forder = as.numeric(poly_smooth_order), dorder = as.numeric(deriv_order))
      # Step required since we lose some wavenumbers due to filter length
      smooth_tmp <- t(smooth_spectra)
      smooth_spectra <- smooth_tmp[,-c(1:((filter_length + 1)/2),
                                       ((ncol(smooth_tmp) - (1+(filter_length + 1)/2)):ncol(smooth_tmp)))]

      spectra <- smooth_spectra
      wavenumber <- data.frame(wavenumber[-c(1:((filter_length + 1)/2),
                                             ((ncol(smooth_tmp) - (1+(filter_length + 1)/2)):ncol(smooth_tmp))),])
    }

    # Normalise selection
    if(norm_meth == "norm_p"){
      spectra <- norm_p(spectra, bin_size = bin_size)
    } else if(norm_meth == "norm_minmax"){
      spectra <- norm_minmax(spectra)
    } else if(norm_meth == "norm_snv"){
      spectra <- norm_snv(spectra)
    } else if(norm_meth == "norm_vec"){
      spectra <- norm_vec(spectra)
    } else if(norm_meth == "norm_a"){
      spectra <- norm_a(spectra, bin_size = bin_size)
    } else if(norm_meth == "none"){
      spectra <- spectra
    }

    return(as.data.frame(spectra))

  }
}
