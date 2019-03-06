#' Standardize the Unrotated Factor Loadings
#'
#' This function standardizes the unrotated factor loadings using two methods: Kaiser's normalization and Cureton-Mulaik standardization.
#'
#' @param method (Character) The method used for standardization. There are three option: "none", "Kaiser", and "CM".
#' \itemize{
#'   \item \strong{"none"}: No standardization is conducted on the unrotated factor loadings matrix
#'   \item \strong{"Kaiser"}: The rows of the unrotated factor loadings matrix are rescaled to have unit-lengths.
#'   \item \strong{"CM"}: Apply the Cureton-Mulaik standardization to the unrotated factor loadings matrix.
#' }
#' @param lambda (Matrix) The unrotated factor loadings matrix (or data frame).
#' 
#' @references Browne, M. W. (2001). An overview of analytic rotation in exploratory factor analysis. \emph{Multivariate Behavioral Research, 36}(1), 111-150.
#' @references Cureton, E. E., & Mulaik, S. A. (1975). The weighted varimax rotation and the promax rotation. \emph{Psychometrika, 40}(2), 183-195.
#' 
#' @family Factor Analysis Routines
#' 
#' @return The resulting output can be used to standardize the factor loadings as well as providing the inverse matrix used to unstandardize the factor loadings after rotating the factor solution.
#' \itemize{
#'   \item \strong{Dv}: (Matrix) A diagonal weight matrix used to standardize the unrotated factor loadings. Pre-multiplying the loadings matrix by the diagonal weight matrix (i.e., Dv %*% lambda) is how to standardization occurs.
#'   \item \strong{DvInv}: (Matrix) The inverse of the diagonal weight matrix used to standardize. To unstandardize the ultimate rotated solution, pre-multiply the rotated factor loadings by the inverse of Dv (i.e., DvInv %*% rotatedLambda) to obtain finish the standardization process.
#'   \item \strong{lambda}: (Matrix) The standardized, unrotated factor loadings matrix.
#'   \item \strong{unstndLambda}: (Matrix) The original, unstandardized, unrotated factor loadings matrix. (DvInv %*% lambda == unstndLambda)
#' }
#' @export

faStandardize <- function(method,
                        lambda) {
  ## Purpose: Find the diagonal weight matrix to standardize lambda
  ##
  ## Args: method: (Character) Which method to use: "none", "Kaiser", "CM".
  ##
  ## Output: Dv: (Matrix) Diagonal weight matrix for standardizing
  ##         DvInv: (Matrix) Inverse of the Dv matrix for later unstandardizing
  ##         lambda: (Matrix) Standardized factor loading matrix
  ##

  #### --------- ERROR CHECKING --------- ####

  ## Check method

  ## Correctly specified?
  if ( method %in% c("none", "Kaiser", "CM") == FALSE ) {
    stop("The 'method' argument is misspecified.")
  } # END ( method %in% c("none", "Kaiser", "CM") == FALSE )

  ## Check lambda

  ## Matrix or data frame?
  if ( class(lambda) %in% c("matrix", "data.frame", "loadings") == FALSE ) {
    stop("The class of 'lambda' must be of class matrix, data.frame, or loadings.")
  } # END if ( class(lambda) %in% c("matrix", "data.frame", "loadings") == FALSE )

  #### --------- DEFINE UTILITY FUNCTIONS --------- ####

  ## Compute Kaiser rotation
  Kaiser <- function(x) {
    diag(sqrt(diag(x %*% t(x)))^-1)
  } # END kaiser

  ## Compute Cureton-Mulaik rotation
  CM <- function(x) {
    NFac <- ncol(x)
    NVar <- nrow(x)
    wghts <- rep(0, NVar)
    fpls <- x[,1] # First principle component loadings
    acosi <- acos( NFac ^ (-1/2) )
    for (i in 1:NVar) {
      num <- (acosi - acos(abs(fpls[i])))
      dem <- (acosi -
                (function(a, m) ifelse(abs(a) < (m ^ (-1/2)),
                                       pi/2, 0) )(fpls[i], NFac))
      wghts[i] <- cos(num / dem * pi / 2) ^ 2
    }
    return( diag(wghts) )
  } # END CM

  #### --------- BODY OF FUNCTION --------- ####

  ## Save unstandardized lambda
  originalLambda <- lambda

  ## Conduct standardization on lambda before rotation
  switch(method,
         "none" = {
           ## If none, this is identity matrix for later unstandardization
           Dv <- DvInv <- diag(nrow(lambda))
         },
         "Kaiser" = {
           ## Diagonal weight matrix to norm rows to have unit length
           Dv <- Kaiser(lambda)

           ## Find inverse for later unstandardization
           DvInv <- solve(Dv)

           ## Standardize (norm) the unrotated loading matrix
           lambda <- Dv %*% lambda
         },
         "CM" = {
           ## Compute Kaiser diagonal weight matrix
           Dk <- Kaiser(lambda)

           ## Compute Cureton-Mulaik D weight matrix on Kaiser-normed lambda
           Dcm <- CM(Dk %*% lambda)

           ## Diagonal matrix to normalize lambda
           Dv <- Dcm %*% Dk

           ## Inverse to later unstandardize rotated solution
           DvInv <- solve(Dv)

           ## Standardize lambda
           lambda <- Dv %*% lambda
         }
  ) # END switch()

  list(Dv           = Dv,
       DvInv        = DvInv,
       lambda       = lambda,
       unstndLambda = originalLambda)

} # END Standardize
