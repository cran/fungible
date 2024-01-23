# Niels Waller
# Updated April 20, 2023
#'
#' Generate the marginal density of a correlation from a uniformly 
#' sampled R matrix.
#' 
#' @param NVar  (integer) The order of the correlation matrix.  
#'
#' @return \code{corDensity} returns the following objects:
#' \itemize{
#'    \item  \strong{r} (numeric) A sequence of numbers from -1, to 1 in .001 
#'    increments.
#'    \item  \strong{rDensity} (numeric) The density of \code{r}.
#' }
#'
#' @references Hürlimann, W.  (2012).  Positive semi-definite correlation 
#' matrices: Recursive algorithmic generation and volume measure.  
#' Pure Mathematical Science, 1(3), 137--149.  
#' 
#' Joe, H.  (2006).  Generating random correlation matrices based on 
#' partial correlations. Journal of Multivariate Analysis, 97(10), 2177--2189.  
#' 
#' @author Niels G. Waller
#'
#' @examples
#'
#'   out <- corDensity(NVar = 5)
#'   
#'   plot(out$r, out$rDensity, 
#'        typ = "l",
#'        xlab = "r",
#'        ylab = "Density of r",
#'        main = "")
#' @export         
#'
corDensity <- function(NVar){
  r <- seq(-1, 1, .001)
  d <- NVar
  dofr <- function(r)
  {
    (1 - r^2)^(d/2 - 1)
  }
  rDensity <- dofr(r)/(.001 * ( sum(dofr(r)) ) )
  list(r = r,
       rDensity = rDensity)
}
