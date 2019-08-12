#' Rolling Circle Filter- Dev - use with apply
#'
#' Background removal for spectral data using RCF, Pchip used for baseline fit
#' @param raw_spec Spectra data that has been interpolated already
#' @param radius RCF radius selection
#' @param wavenumber Wavenumber range included
#' @return Corrected spectra with RCF applied (background removed)
#' @export
#' @import pracma

RCF_GENERALISED <- function(raw_spectra, radius, wavenumber){

  rad <- as.numeric(radius)

  # Creating our semicircle for RCF Interval is automatically set to 1
  x <- seq(-rad, rad)
  y <- sqrt(rad^2 - x^2)

  ## Creating artificial start point for the circle
  # X coordinates

  # Creating our new x data for the circle to roll under
  fakexstart <- wavenumber[1,] -  seq(from = 1, to = length(x))
  fakexend <- wavenumber[nrow(wavenumber),] + seq(from = 1, to = length(x))
  # And now the entire set of x data - sandwiches the actual wavenumber range
  # between the fake start data, and fake end data (for x axis)
  totalxdata <- c(rev(fakexstart), wavenumber[,1], fakexend)

  # Y coordinates
  # Looks for the max y in each spectral set so we don't remove them!
  fakeycoords <- max(raw_spectra)

    # again sandwiching actual y data between the fake y data
  totalycoords <- rep(fakeycoords, length(fakexstart))
  totalydata <- c(t(totalycoords), raw_spectra, t(totalycoords))
  # Now combinding the x and y test data sets together for each spectra
  testdataset <- rbind(totalxdata, totalydata)

  ### Aspect ration calculation
  ## So the circle scales to the size of the data
  # and adjusts start points for x and y
  AspectRatio <- 2*(max(raw_spectra)/(max(wavenumber)-min(wavenumber)))
  AspectYCoords <- AspectRatio*y
  CircleYStartpoint <- AspectYCoords + as.numeric(fakeycoords) - AspectYCoords[[ceiling((length(x)/2))]]

  CirclestartpointX <- as.data.frame(x+totalxdata[1]+ rad)

  CircleStartY <- cbind(as.data.frame(CirclestartpointX), as.data.frame(CircleYStartpoint))

  ## Data in circle


  # Loops over all spectra and runs the RCF on each in turn
  loopsize <- seq(from = 1, to = (ncol(testdataset) - length(x)), by = 1)
  b <- length(x)
  indices <- loopsize-1+b

  idx <- seq2(from = loopsize,
       to = indices, by = 1)
  goal <- list()
  for (i in 1:ncol(idx)){
    goal[[i]] <- testdataset[, idx[,i]]
  }
  Sys.time() - t1

  dist <- matrix(ncol = length(goal), nrow = b)
  for (i in 1:length(goal)){
    dist[,i] <- as.numeric(goal[[i]][2,]) - CircleStartY[,2]
  }
  distance <- t(dist)

  Min <- matrix(ncol = 2, nrow = nrow(distance))
  for (i in 1:nrow(distance)){
    Min[i,] <- c(match(min(distance[i,]), distance[i,]), min(distance[i,]))
  }

  xbit <- matrix(nrow = nrow(Min), ncol = 2)
  for (i in 1:nrow(Min)){
    xbit[i,] <- unlist(goal[[i]][,Min[i,1]])
  }
  XY <- xbit
  Unique <- unique(XY)
  Xaxis <- Unique[,1]
  Yaxis <- Unique[,2]

  baseline <- pchip(Xaxis,Yaxis,wavenumber[,1])
  corrected <- t(as.data.frame(raw_spectra)) - as.data.frame(baseline)


  corrected[corrected<0] <- 0
  return(t(corrected))

}

