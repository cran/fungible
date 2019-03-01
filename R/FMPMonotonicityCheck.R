## January 27, 2016
## updated February 15, 2016
## This function is used to determine whether the derivative of the 
## polynomial is everywhere positive. in FMP models the derivative of the polynomial
## must be positive everywhere to insure monotonicity.  


#' Utility function for checking FMP monotonicity
#' 
#' Utility function for checking whether candidate FMP coefficients yield a
#' monotonically increasing polynomial.
#' 
#' 
#' @param b A vector of 8 polynomial coefficients (\eqn{b}) for
#' \eqn{m(\theta)=b_0 + b_1 \theta + b_2 \theta^2 + b_3 \theta^3 + b_4 \theta^4
#' + b_5 \theta^5 + b_6 \theta^6 + b_7 \theta^7}.
#' @param lower,upper \eqn{\theta} bounds for monotonicity check.
#' @param PLOT Logical (default = FALSE).  If PLOT = TRUE the function will
#' plot the original polynomial function for \eqn{\theta} between lower and
#' upper.
#' @return \item{increasing}{Logical indicating whether function is
#' monotonically increasing.} \item{minDeriv}{Minimum value of the derivative
#' for the polynomial.} \item{minTheta}{Value of \eqn{\theta} at derivative
#' minimum.}
#' @author Niels Waller
#' @keywords statistics
#' @importFrom graphics plot
#' @export
#' @examples
#' 
#' 
#' ## A set of candidate coefficients for an FMP model.
#' ## These coefficients fail the test and thus
#' ## should not be used with genFMPdata to generate
#' ## item response data that are consistent with an 
#' ## FMP model.
#'  b <- c(1.21, 1.87, -1.02, 0.18, 0.18, 0, 0, 0)
#'  FMPMonotonicityCheck(b)
#' 
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
