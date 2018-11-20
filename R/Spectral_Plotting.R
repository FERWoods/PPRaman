#' Spectral Plotting
#'
#' Easy and clear plotting for spectral data for processed and pre-processed spectra
#' @param raw_spec Spectra in matrix form (like raw spectra)
#' @param samp_num # of samples
#' @return plot of spectra
#' @export
#'

spectral_plots <- function(raw_spec, samp_num, wavenumber){

  OBJ_For_Plotting <- t(raw_spec)
  #Add the wavenumber to the data
  OBJ_For_Plotting <- cbind(wavenumber, OBJ_For_Plotting)

  #change colnames
  samples <- seq(from = 1, to = samp_num, by = 1)
  samplesno <- paste("Sample", samples, sep = "" )

  colnames(OBJ_For_Plotting) <- c("Wavenumber", samplesno)

  # use melt function to make plotting easier to include different samples in 1 plot
  OBJ_For_Plotting <- melt(OBJ_For_Plotting, id.vars = "Wavenumber")

  # Creates basis of plot
  ptitle <- readline("Title for Spectral Plot?")
  p <- ggplot(OBJ_For_Plotting) + geom_line(aes(Wavenumber,value, col = variable)) +
    xlab(expression(paste("Wavenumber (", cm^-1, ")" ))) + ylab("Intensity (a.u)")
    ggtitle(ptitle) + theme(plot.title = element_text(hjust = 0.5))

  return(p)
}
