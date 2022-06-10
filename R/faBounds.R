#' Bounds on the Correlation Between an External Variable and a Common Factor
#' 
#'
#' This function computes the  bounds on the correlation between an 
#' external variable and a common factor.
#'
#' @param  Lambda  (matrix) A p x 1 matrix of factor loadings. 
#' @param  RX      (matrix) A p x p matrix of correlations for the factor indicators.
#' @param  rXY     (vector) A p x 1 vector of correlations between the factor
#' indicators (X) and the external variable (Y).
#' @param  alphaY  (scalar)  The reliability of Y. Default \code{alphaY = 1}.
#'
#' @return \code{faBounds} returns the following objects:
#' \itemize{
#'   \item \strong{Lambda} (matrix) A p x 1 vector of factor loadings. 
#'    \item  \strong{RX} (matrix) The indicator correlation matrix. 
#'    \item  \strong{rXY}: (vector) The correlations between the factor indicators (X) and the 
#'    external variable (Y).
#'    \item \strong{alphaY} (integer) The reliability of the external variable.
#'    \item \strong{bounds} (vector)  A 2 x 1 vector that includes the lower and upper bounds 
#'    for the correlation between an external variable and a common factor. 
#'    \item \strong{rUiY} (vector) Correlations between the unique factors and the 
#'    external variable for the lower bound estimate.
#'    \item \strong{rUjY} (vector) Correlations between the unique factors and the 
#'    external variable for the upper bound estimate.
#' }
#'
#' @references Steiger, J. H.  (1979).  The relationship between external 
#' variables and common factors. Psychometrika, 44, 93-97. 
#' 
#' @references Waller, N. G. (under review). New results on the relationship 
#' between an external variable and a common factor. 
#' 
#' @author Niels G. Waller
#'
#' @examples
#' ## Example 
#' ## We wish to compute the bounds between the Speed factor from the 
#' ## Holzinger (H) and Swineford data and a hypothetical external 
#' ## variable, Y.
#' 
#' ## RH = R matrix for *H*olzinger Swineford data
#' RH <- 
#'  matrix(c( 1.00,   0,    0,     0,     0,     0,
#'            .73, 1.00,    0,     0,     0,     0, 
#'            .70,  .72,  1.00,    0,     0,     0,
#'            .17,  .10,   .12,  1.00,    0,     0,
#'            .11,  .14,   .15,   .49,  1.00,    0,
#'            .21,  .23,   .21,   .34,   .45,  1.00), 6, 6)
#'
#' RH <- RH + t(RH) - diag(6)
#' RX <- RH[4:6, 4:6]
#'
#' ## S-C = Straight-curved
#'  colnames(RX) <- rownames(RX) <-
#'         c("Addition", "Counting dots", "S-C capitals")
#' print( RX, digits = 2 ) 
#'
#' ## Extract 1 MLE factor  
#' fout <- faMain(R = RX, 
#'               numFactors = 1, 
#'               facMethod = "faml", 
#'               rotate="none")
#'
#' ## Lambda = factor loadings matrix  
#' Lambda <- fout$loadings
#' print( Lambda, digits = 3 ) 
#' 
#' ## rXY = correlations between the factor indicators (X) and
#' ## the external variable (Y)
#'
#'  rXY = c(.1, .2, .3)
#'  
#'  # Assume that the reliability of Y = .75
#'  
#'  faBounds(Lambda, RX, rXY, alphaY = .75)
#'  
#' @export
faBounds <- function(Lambda, RX, rXY, alphaY = 1){
  
  # h2 = commonalities of factor indicators
  h2 <- diag(Lambda %*% t(Lambda))
  
  # Correlations between X and the Unique factors
  Psi <- diag(sqrt( 1 - h2 ))
  
  # Number of factor indicators
  p <- ncol(RX)
  
  Psi.inv <- solve(Psi)
  RXinv <- solve(RX)
  
  # Disattenuate rXY  for Y unreliability
  rXY <- rXY/sqrt(alphaY)
  
  # cov_F_Y.X = the covariance of F and Y.X (Y.X = the  
  # part of Y that is orthogonal to X)
  cov_F_Y.X <- t(Lambda) %*% RXinv %*% rXY
  
  # sig1 = std dev of Y.X
  sig1 <- (1 - t(rXY) %*% RXinv %*% rXY)^.5
  
  # sig2 = std dev F.X
  sig2 <- (1 - t(Lambda) %*% RXinv %*% Lambda)^.5
  
  # Compute bounds taking the unreliability of Y into account
  LB <- (cov_F_Y.X - sig1 * sig2) * sqrt(alphaY)
  UB <- (cov_F_Y.X + sig1 * sig2) * sqrt(alphaY)
  
  # attenuate rXY to original values
  rXY <- rXY * sqrt(alphaY)
  
  # Compute correlations between the unique factors and the external variable
  rUiY <- Psi.inv %*% (rXY - Lambda * as.numeric(LB))
  rUjY <- Psi.inv %*% (rXY - Lambda * as.numeric(UB))
  
  
  ## RETURN
  rownames(rUiY) <- rownames(rUjY) <-paste0("U", 1:p)
  colnames(rUiY) <- colnames(rUjY) <- "Y"
  colnames(Lambda) <- "F"
  rownames(Lambda) <- paste0("X", 1:p)
  

  list( Lambda = Lambda,
        RX = RX,
        rXY = rXY,
        alphaY = alphaY,
        bounds = c(LB,UB),
        rUiY   = rUiY,  # rUY for LB
        rUjY   = rUjY)  # rUY for UB
} #END faBounds
