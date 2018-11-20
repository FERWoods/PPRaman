#' Rolling Circle Filter
#'
#' Background removal for spectral data using RCF
#' @param raw_spec Spectra data that has been interpolated already
#' @param samp_num Number of samples within data - user input
#' @param baseline_fit Choice of baseline fit from Pchip or Linear
#' @param radius RCF radius selection
#' @return Corrected spectra with RCF applied (background removed)
#' @export

RCF <- function(raw_spec, samp_num, baseline_fit, radius){

  raw_spectra <- raw_spec
  samp_num <- as.numeric(samp_num)
  rad <- radius
  # Creating our semicircle for RCF Interval is automatically set to 1
  x <- seq(-rad, rad)
  y <- sqrt(rad^2 - x^2)

  ## Creating artificial start point for the circle
  #X coordinates
  fakexstart <- list()
  fakexend <- list()
  for (i in 1:length(x)){
    fakexstart[i] <- wavenumber[1,] - i
    fakexend[i] <- wavenumber[nrow(wavenumber),] + i
  }
  # And now the entire set of x data - sandwiches the actual wavenumber range
  # between the fake start data, and fake end data (for x axis)
  totalxdata<- as.data.frame((c(rev(fakexstart), as.list(wavenumber[,1]), fakexend)))

  #Y coordinates
  # this just looks for the max y in each spectral set so we don't remove them!
  fakeycoords <- list()
  totalycoords <- list()
  for (i in 1:samp_num){
    fakeycoords[i] <- max(raw_spectra[i,])
  }

  # again sandwiching actual y data between the fake y data
  totalycoords <- as.data.frame(replicate(n = length(fakexstart), fakeycoords))
  totalydata <- as.data.frame(cbind(totalycoords, raw_spectra, totalycoords))

  # Now combinding the x and y test data sets together for each spectra
  testdataset <- list()
  colnames(totalxdata) <- colnames(totalydata)
  for (i in 1:samp_num){
    testdataset[[i]] <- rbind(totalxdata, totalydata[i,])
  }


  ### Aspect ration calc
  ## So the circle scales to the size of the data
  # and adjusts start points for x and y
  AspectRatio <- matrix(ncol = samp_num, nrow = 1)
  AspectYCoords <- matrix(ncol = length(x), nrow = samp_num)
  CircleYStartpoint <- matrix(ncol = length(x), nrow = samp_num)
  for (i in 1:samp_num){
    AspectRatio[,i] <- 2*(max(raw_spectra[i,])/(max(wavenumber)-min(wavenumber)))
    AspectYCoords[i,] <- AspectRatio[,i]*y
    CircleYStartpoint[i,] <- AspectYCoords[i,] + as.numeric(fakeycoords[i]) - AspectYCoords[i,ceiling((length(x)/2))]
  }

  CirclestartpointX <- as.data.frame(x+totalxdata[1,1]+ rad)

  CircleStartY <- list()
  for (i in 1:samp_num){
    Circlestartpos <- cbind(as.data.frame(CirclestartpointX), as.data.frame(CircleYStartpoint[i,]))
    CircleStartY[[i]] <- Circlestartpos
  }

  ## Data in circle
  corrected <- list()

  # Loops over all spectra and runs the RCF on each in turn
  for (j in 1:samp_num){
    loopsize <- ncol(testdataset[[1]]) - length(x)
    b <- length(x)

    goal <- list()
    for (i in 1:loopsize){
      ind <- i-1+b
      goal[[i]] <- as.data.frame(t(testdataset[[j]][,i:ind]))
    }

    dist <- matrix(ncol = length(goal), nrow = b)
    for (i in 1:length(goal)){
      dist[,i] <- as.numeric(goal[[i]][,2]) - CircleStartY[[j]][,2]
    }
    distance <- t(dist)

    Min <- matrix(ncol = 2, nrow = nrow(distance))
    for (i in 1:nrow(distance)){
      Min[i,] <- c(match(min(distance[i,]), distance[i,]), min(distance[i,]))
    }

    xbit <- matrix(nrow = nrow(Min), ncol = 2)
    for (i in 1:nrow(Min)){
      xbit[i,] <- unlist(goal[[i]][Min[i,1],])
    }
    XY <- xbit
    Unique <- unique(XY)
    Xaxis <- Unique[,1]
    Yaxis <- Unique[,2]

    baseline <- pchip(Xaxis,Yaxis,wavenumber[,1])

    corrected[[j]] <- t(as.data.frame(raw_spectra[j,])) - as.data.frame(baseline)


  }
  corrected_spec <- do.call("cbind", corrected)
  corrected_spec[corrected_spec<0] <- 0
  return(corrected_spec)

}
