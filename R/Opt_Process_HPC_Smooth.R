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
#' @return Preprocessed spectra
#' @export


opt_process_hpc_smooth <- function(spectra, norm_meth, rm_bl, RCF_rad, poly_bl_order,
                                filter_length, deriv_order, poly_smooth_order, smoothing_opt){
  if(smoothing_opt == 1){
    smooth_spectra <- apply(spectra, 1, savgol, fl = as.numeric(filter_length),
                            forder = as.numeric(poly_smooth_order), dorder = as.numeric(deriv_order))
    # Step required since we lose some wavenumbers due to filter length
    smooth_spectra <- t(smooth_spectra)
    smooth_spectra <- smooth_spectra[,-(1:((filter_length + 1)/2))]
    smooth_spectra <- smooth_spectra[,-((ncol(smooth_spectra) - (1+(filter_length + 1)/2)):ncol(smooth_spectra))]
    spectra <- smooth_spectra
    wavenumber <- data.frame(wavenumber[-c((1:((filter_length + 1)/2)),
                                           ((ncol(smooth_spectra) - (1+(filter_length + 1)/2)):ncol(smooth_spectra))),])
  }

  # Baseline removal selection
  if(rm_bl == "rcf"){
    bl_rmv <- apply(spectra, 1, RCF_GENERALISED, radius = as.numeric(RCF_rad), wavenumber = wavenumber)
  } else if(rm_bl == "der"){
    bl_rmv <- spectra
  } else if(rm_bl == "rub"){
    temp <- spc.rubberband(spectra)
    bl_rmv <- spectra - temp
  } else if(rm_bl == "pol"){
    bl_rmv <- poly_baseline_removal(as.hyperSpec(spectra), poly.order = poly_bl_order)
    #smooth_replacement <- zeros(m = (3 + 1)/2, n = nrow(bl_rmv))
    #smooth_replacement2 <- cbind(smooth_replacement, spectra, smooth_replacement)
    #bl_rmv <- smooth_replacement2
  } else if(rm_bl == "none"){
    bl_rmv <- spectra
  }

  # Normalise selection
  if(norm_meth == "norm_p"){
    norm_spec <- norm_p(t(bl_rmv))
  } else if(norm_meth == "norm_minmax"){
    norm_spec <- norm_minmax(t(bl_rmv))
  } else if(norm_meth == "norm_snv"){
    norm_spec <- norm_snv(t(bl_rmv))
  } else if(norm_meth == "norm_vec"){
    norm_spec <- norm_vec(t(bl_rmv))
  } else if(norm_meth == "norm_a"){
    norm_spec <- norm_a(t(bl_rmv))
  } else if(norm_meth == "none"){
    norm_spec <- t(bl_rmv)
  }

  return(as.data.frame(norm_spec))
}


