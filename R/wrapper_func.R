#' Options for pre-processing spectra - developing version
#'
#' @param spectra Spectra data that has been interpolated already
#' @param norm_meth Choice of normalisation
#' @param rm_bl Choice of background removal
#' @param baseline_fit Baseline fit for RCF - currently off
#' @param RCF_rad RCF radius selection
#' @param poly_order Choice of polynomial order for SavGol etc
#' @return Preprocessed spectra
#' @export

opt_process_hpc_dev <- function(spectra, norm_meth, rm_bl, baseline_fit, RCF_rad, poly_order,
                                  filter_length, deriv_order){

  # Baseline removal selection
  if(rm_bl == "rcf"){
      bl_rmv <- apply(spectra, 1, RCF_dev, baseline_fit, as.numeric(RCF_rad))
  } else if(rm_bl == "der"){
      bl_rmv <- apply(spectra, 1, savgol, fl = as.numeric(filter_length),
                      forder = as.numeric(poly_order), dorder = as.numeric(deriv_order))
      # Step required since we lose some wavenumbers due to filter lenght
      bl_rmv[1:(filter_length + 1)/2, ] <- 0
      bl_rmv[(ncol(spectra) - (filter_length + 1)/2):ncol(spectra), ] <- 0
  } else if(rm_bl == "rub"){
      temp <- spc.rubberband(spectra)
      bl_rmv <- spectra - temp
  } else if(rm_bl == "pol"){
      bl_rmv <- poly_baseline_removal(spectra, poly.order = poly_order)
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

