#' Reading in raw spectra and interpolating
#'
#' Reads in mutiple raw spectra and interpolates using hyperSpec pkg
#' @return Interpolated spectra in hyperSpec format
#' @export

read_interp_spectra <- function(){
  list.files()
  filelist = list.files(pattern = ".*.txt")

  #assuming tab separated values with a header
  datalist = lapply(filelist, function(x)read.table(x, header=T))

  #assuming the same header/columns for all files
  datafr = do.call("rbind", datalist)
  files <- choose.files()
  # This reads all the files in using read.table
  inputfiles <- lapply(files, read.table)

  #Storing the original wavenumber for creating hyperspec object
  old_wn <- inputfiles[[1]][,1]

  # remove the w/n column from each spectra
  rmv_wn <- lapply(inputfiles, function(x) { x["V1"] <- NULL; x })


  raw_spec <- t(do.call("cbind", rmv_wn))

  #Step sizes for interpolating
  waveno <- seq(611.6, 1717, by = 1.09)


  #Converts raw spectra to hyperspec object to work with the package
  raw_hyperSpec <- as.hyperSpec(raw_spec, wl = old_wn)

  #Add labels

  raw_hyperSpec@label$.wavelength <- expression(Wavenumber (cm^-1))
  raw_hyperSpec@label$spc <- expression(Intensity (a.u))

  #Interpolating the data using loess
  int_spec <- spc.loess(raw_hyperSpec, waveno)

  return(int_spec)
}


