# Not run / tested -- this is another method for shifting the spectra
# based on Claudia B's paper -- similar method
# uses square polynomial, received error for 4th order, same issue with polyfit func.
find.max <- function (y, x){
  pos <- which.max (y) + (-1:1)
  X <- x [pos] - x [pos [2]]
  Y <- y [pos] - y [pos [2]]

  X <- cbind (1, X, X^2) # fits to a square polynomial
  coef <- qr.solve (X, Y)

  - coef [2] / coef [3] / 2 + x [pos [2]]
}



#bandpos <- apply (Blood_Spectra[[,, 800 ~ 880]], 1, find.max, wl (Blood_Spectra [,, 800 ~ 880]))
#refpos <- find.max (colMeans (Blood_Spectra[[,, 800 ~ 880]]), wl (Blood_Spectra [,, 800 ~ 880]))
#shift1 <- refpos - bandpos
