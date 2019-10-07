#' Ledermann's inequality for factor solution identification
#'
#' Ledermann's (1937) inequality to determine either (a) how many factor 
#' indicators are needed to uniquely estimate a user-specified number 
#' of factors or (b) how many factors can be uniquely estimated from 
#' a user-specified number of factor indicators. See the \strong{Details} 
#' section for more information
#'
#' @param numFactors (Numeric) Determine the number of variables needed
#' to uniquely estimate the [user-specifed] number of factors. Defaults 
#' to \code{numFactors = NULL}. 
#' @param numVariables (Numeric) Determine the number of factors that can be
#' uniquely estimated from the [user-specifed] number of variables Defaults 
#' to \code{numVariables = NULL}. 
#' 
#' @details The user will specified either (a) \code{numFactors} or (b) 
#' \code{numVariables}. When one value is specified, the obtained estimate 
#' for the other may be a non-whole number. If estimating the number of 
#' required variables, the obtained estimate is rounded up 
#' (using \code{\link[base]{ceiling}}). If estimating the number of factors,
#' the obtained estimate is rounded down (using \code{\link[base]{floor}}). For example,
#' if \code{numFactors = 2}, roughly 4.56 variables are required for an identified
#' solution. However, the function returns an estimate of 5.
#' 
#' @details For the relevant equations, see Thurstone (1947, p. 293) Equations 10 
#' and 11.  
#' 
#' @return 
#' \itemize{
#'    \item \bold{numFactors} (Numeric) Given the inputs, the number of factors 
#'    to be estimated from the \code{numVariables} number of factor indicators. 
#'    \item \bold{numVariables} (Numeric) Given the inputs, the number of 
#'    variables needed to estimate \code{numFactorso}. 
#' }
#' 
#' @references Ledermann, W. (1937). On the rank of the reduced correlational 
#' matrix in multiple-factor analysis. \emph{Psychometrika, 2}(2), 85-93.
#' @references Thurstone, L. L. (1947). Multiple-factor analysis; a development and expansion of The Vectors of Mind.
#' 
#' @family Factor Analysis Routines
#' 
#' @author Casey Giordano
#' 
#' @examples 
#' ## To estimate 3 factors, how many variables are needed?
#' Ledermann(numFactors   = 3,
#'           numVariables = NULL) 
#'           
#' ## Provided 10 variables are collected, how many factors 
#'   ## can be estimated?
#' Ledermann(numFactors   = NULL,
#'           numVariables = 10)
#' 
#' @export

Ledermann <- function(numFactors = NULL, 
                      numVariables = NULL) {

  ## References: Ledermann (1937)
  ## Thurstone (1947, p. 293) equations 10 & 11
  if (is.null(numFactors) && is.null(numVariables)) {
    stop("Must specify either 'numFactors' or 'numVariables'.")
  } # END if (numFactors == NULL && numVariables == NULL)


  if (is.null(numFactors) && is.numeric(numVariables)) {
    ## Thurstone (1947, p. 293) equation 10
    numFactors <- ((2 * numVariables + 1) - sqrt(8 * numVariables + 1)) / 2
    
    
    
    if(floor(numFactors) == 1){
      cat("The maximum number of factors to be estimated from",
          numVariables, "variables is", floor(numFactors), "factor.\n\n")
    }
    else{
    cat("The maximum number of factors to be estimated from",
        numVariables, "variables is", floor(numFactors), "factors.\n\n")
    }
    
    
  } # END if (numFactors == NULL && is.numeric(numVariables))

  if (is.null(numVariables) && is.numeric(numFactors)) {
    ## Thurstone (1947, p. 293) equation 11
    numVariables <- ((2 * numFactors + 1) + sqrt(8 * numFactors + 1)) / 2

    cat("The minimum number of variables needed to identify",
        numFactors, "factors is", ceiling(numVariables), "variables.\n\n")
  } # END if (numVariables == NULL && is.numeric(numFactors))

  list(numFactors   = floor(numFactors),
       numVariables = ceiling(numVariables))
} # END Ledermann <- function(numFactors = NULL, numVariables = NULL)

