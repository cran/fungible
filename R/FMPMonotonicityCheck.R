## January 27, 2016
## updated February 15, 2016
## This function is used to determine whether the derivative of the 
## polynomial is everywhere positive. in FMP models the derivative of the polynomial
## must be positive everywhere to insure monotonicity.  
FMPMonotonicityCheck<-function(b, lower = -20, upper = 20, PLOT = FALSE){
  if(length(b)!=8) stop("\n\nSerious Error: b should have 8 elements\n\n")
  ## b is an 8x1 vector of FMP polynomial coefficients
  ## b0, b1, . . . b7
  DerivFMP <-  function(theta,b) {
      b[2] + 2*b[3]*theta + 3*b[4]*theta^2 + 4*b[5]*theta^3 + 5*b[6]*theta^4 +
      6*b[7]*theta^5 + 7*b[8]*theta^6
  }
  xx<-optimize(f = DerivFMP, interval = c(lower,upper), b)
  
  if( is.na(xx$objective) ) stop("\n\n\nSerious optimization problems occured\n\n")
  
  
  cat(paste("\nMinimum derivative of the polynomial in the interval (", lower, ",", upper,") = ", 
            round(xx$objective,4), sep=""))
  cat("\nThe minimum derivative occurs at theta = ", round(xx$minimum,4))
  
 
  
  monotonic <- ( xx$objective > 0 )
  if( !monotonic) cat("\nWarning: Polynomial is not monotonically increasing.\n\n")
 
  minDeriv<-xx$objective
  
  if(PLOT){
    x <- seq(lower,upper, by=.01)
    X <- cbind(1,x,x^2, x^3, x^4, x^5, x^6, x^7)
    y <- X %*% b
    plot(x, y, type="l", lwd=2, col="blue", 
         xlab=expression(theta),
         ylab=expression(f(theta)))
  }
  
  list(increasing = monotonic, minDeriv = minDeriv, minTheta = xx$minimum)	  
  }
