# Niels Waller
# Updated April 19, 2023
#'
#' Compute the volume of the elliptope of possible correlation 
#' matrices of a given dimension.
#' 
#'
#'
#' @param NVar  (integer) The size of each correlation matrix in the elliptope. 
#' For instance, if we are interested in the volume of the space of all possible 
#' 5 x 5  correlation matrices then \code{NVar = 5}. 
#'
#' @return \code{VolElliptope} returns the following objects:
#' \itemize{
#'    \item  \strong{VolElliptope} (numeric) The volume of the elliptope.
#'    \item  \strong{VolCube}: (numeric) The volume of the embedding 
#'    hyper-cube.
#'    \item \strong{PrcntCube} (numeric) The percent of the hyper-cube that is 
#'    occupied by the elliptope. \code{PrcntCube = 100 x VolElliptope/VolCube}.
#' }
#'
#' @references Joe, H. (2006). Generating random correlation matrices 
#' based on partial correlations. *Journal of Multivariate Analysis*, 
#' *97* (10), 2177--2189. 
#' 
#' Hürlimann, W. (2012). Positive semi-definite 
#' correlation matrices: Recursive algorithmic generation and volume 
#' measure. *Pure Mathematical Science, 1* (3), 137--149. 
#' 
#' @author Niels G. Waller
#'
#' @examples
#' # Compute the volume of a 5 x 5 correlation matrix.
#' 
#' VolElliptope(NVar = 5)
#'
#' @export 
VolElliptope <- function(NVar){
  
  # Compute the volume of the elliptope of PSD R matrices
  
  VolE <- 1
  
  for(k in 1:(NVar - 1)){
    f <- function(t) (1 - t^2)^(.5*(k-1))  
    Jk <-   integrate(f,
                      lower = -1, 
                      upper =  1)$value
    VolE <- VolE * Jk^k
  }
  
  VolCube <- (2)^(.5 * NVar * (NVar - 1))
  PrcntCube = 100 * VolE/VolCube
  
  list(VolElliptope = VolE,
       VolCube = VolCube,
       PrcntCube = PrcntCube)
}

