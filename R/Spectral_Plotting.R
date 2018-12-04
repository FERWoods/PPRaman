#' Spectral Plotting
#'
#' Easy and clear plotting for spectral data for processed and pre-processed spectra
#' @param raw_spec Spectra in matrix form (like raw spectra)
#' @param labels Labels for data
#' @return plot of spectra
#' @export
#' @import ggplot2
#' @import reshape2

spectral_plots <- function(raw_spec, labels){

  OBJ_For_Plotting <- t(raw_spec)

  #Add the wavenumber to the data
  OBJ_For_Plotting <- cbind(wavenumber, OBJ_For_Plotting)

  #change colnames
  label <- labels

  colnames(OBJ_For_Plotting) <- c("Wavenumber", label)

  # use melt function to make plotting easier to include different samples in 1 plot
  OBJ_For_Plotting <- melt(OBJ_For_Plotting, id.vars = "Wavenumber")

  # Creates basis of plot
  #ptitle <- readline("Title for Spectral Plot?")
  p <- ggplot(OBJ_For_Plotting) + geom_line(aes(Wavenumber,value, col = variable)) +
    xlab(expression(paste("Wavenumber (", cm^-1, ")" ))) + ylab("Intensity (a.u)") +
    labs(color = "Sample")
    #ggtitle(ptitle) + theme(plot.title = element_text(hjust = 0.5))

  return(p)
}
