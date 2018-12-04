#' Old version of preprocess seq - smooths too, but not sure it shifts


process_all_spec_historic <- function(){

  # Read in all spectra -- list all in one file
  spectral_info <- read_interp_spectra()
  spectra <- spectral_info[[1]]
  supplementary <- spectral_info[[2]]

  # Run RCF on data
  bl_rmv <- RCF(spectra, nrow(spectra) , "pchip", 150)

  # Normalise to phenylalanine peak
  norm_spec <- norm_p(t(bl_rmv))

  return(list(norm_spec, supplementary))
}
