#' Options V1 Norm & RF Radius

#' First run of options for the HPC -- testing phase
#' @export
#'

# Read in all spectra -- list all in one file
read_in <- function(){
  spectral_info <- read_s
  spectra <- spectral_info[[1]]
  supplementary <- spectral_info[[2]]

  return(list(spectra, supplementary))
}

opt_process_hpc <- function(spectra, norm_meth, RCF_rad){

  # Run RCF on data
  bl_rmv <- do.call(cbind, apply(raw_spec, 1, RCF_dev, 6,
                                   "pchip", 150))


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



