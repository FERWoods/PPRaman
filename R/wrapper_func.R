#' Options for pre-processing spectra - developing version
#'
#' @param spectra Spectra data that has been interpolated already
#' @param norm_meth Choice of normalisation
#' @param RCF_rad RCF radius selection
#' @return Preprocessed spectra
#' @export

#setwd("p:/documents/PPRaman")
opt_process_hpc_dev <- function(spectra, norm_meth, rm_bl, samp_num, baseline_fit, poly_order){

  # Baseline removal selection
  if(rm_bl == "rcf"){
    bl_rmv <- RCF(spectra, as.numeric(samp_num), baseline_fit, as.numeric(poly_order))
  } else if(rm_bl == "der"){
      bl_rm <- SavGol(spectra, poly_order = poly_order)
  } else if(rm_bl == "rub"){
      temp <- spc.rubberband(spectra)
      bl_rm <- spectra - temp
  } else if(rm_bl == "pol"){
      bl_rm <- poly_baseline_removal(spectra, poly.order = poly_order)
  } else if(rm_bl == "none"){
      bl_rm <- spectra
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
  } else if(norm_meth == "none"){
    norm_spec <- t(bl_rmv)
  }

  return(norm_spec)
}
