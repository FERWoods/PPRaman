#' Options for pre-processing spectra - DEV
#'
#' @param spectra Spectra data that has been interpolated already
#' @param norm_meth Choice of normalisation
#' @param RCF_rad RCF radius selection
#' @return Preprocessed spectra
#' @export


opt_process_hpc_dev <- function(spectra, norm_meth, RCF_rad){

  # Run RCF on data
  bl_rmv <- do.call(cbind, apply(spectra, 1, RCF_dev,
                                   "pchip", as.numeric(RCF_rad)))


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
