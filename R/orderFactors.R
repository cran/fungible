

#' Order factor-loadings matrix by the sum of squared factor loadings
#' 
#' Order the columns of a factor loadings matrix in descending order based on 
#' the sum of squared factor loadings. 
#'
#' @param Lambda (Matrix) Factor loadings matrix to be reordered. 
#' @param PhiMat (Matrix, NULL) Factor correlation matrix to be reordered. 
#' @param salient (Numeric) Indicators with loadings < \code{salient} will be 
#' suppressed when computing the factor sum of squares values. Defaults to 
#' salient = .29.
#' @param reflect (Logical) If true, negatively-keyed factors will be reflected.
#' Defaults to reflect = TRUE.
#'
#' @return Returns the sorted factor loading and factor correlation matrices. 
#' \itemize{
#'   \item \strong{Lambda}: (Matrix) The sorted factor loadings matrix.
#'   \item \strong{Phi}: (Matrix) The sorted factor correlation matrix.
#' }
#' 
#' @family Factor Analysis Routines
#' 
#' @examples 
#' \dontrun{
#' Loadings <- 
#'   matrix(c(.49, .41, .00, .00,
#'            .73, .45, .00, .00,
#'            .47, .53, .00, .00,
#'            .54, .00, .66, .00,
#'            .60, .00, .38, .00,
#'            .55, .00, .66, .00,
#'            .39, .00, .00, .68,
#'            .71, .00, .00, .56,
#'            .63, .00, .00, .55), 
#'          nrow = 9, ncol = 4, byrow = TRUE)
#'          
#' fungible::orderFactors(Lambda = Loadings,
#'                         PhiMat = NULL)$Lambda
#' }
#' 
#' @export

orderFactors <- function(Lambda,
                         PhiMat,
                         salient = .29,
                         reflect = TRUE) {
  ## Purpose: Sort the factors into conventional order
  ##
  ## Args: Lambda: (Matrix) Factor loading matrix
  ##       PhiMat: (Matrix) Factor correlation matrix
  
  ## Wipe out low loadings (based on salient argument)
  F0 <- Lambda
  F0[abs(Lambda) < salient] <- 0
  
  ## If only 1 factor, only need to reflect the factor
  if ( ncol(Lambda) == 1 ) {

    ## If factor is negative, reverse it and end the function
    if ( sum(F0) < 0) Lambda <- Lambda * -1
    return( 
      list(Lambda = Lambda,
           PhiMat = PhiMat)
    )
  } # END if ( ncol(Lambda) == 1 )   
  
  ## If no Phi is specified, carry this through to the end
  Phi <- NULL
  
  ## Flip poorly keyed factors (and associated factor correlations)
  if (reflect == TRUE) {
    ## Add the columns, detect the sign. Used to flip negative factors
    FSignMat <- diag(sign(colSums(F0)))
    
    ## Flip negatively keyed factors
    Lambda   <- Lambda %*% FSignMat
    
    if ( !is.null(PhiMat) ) {
      
      ## Also flip Phi matrix if applicable
      Phi      <- FSignMat %*% PhiMat %*% FSignMat
      
    } # END if ( !is.null(PhiMat) ) 
    
  } # END if (reflect == TRUE) 
  
  ## Vector used to re-order the factor loading/correlation matrices
  ## Based on sum of squared loadings
  LambdaSS <- sort.list( colSums(F0^2), decreasing = TRUE )
  
  ## Sort the keyed factor loadings
  Lambda <- Lambda[, LambdaSS]
  
  ## If Phi is specified, sort it
  if ( !is.null(PhiMat) ) {
    
    ## Sort the keyed factor correlations
    I   <- diag(ncol(Lambda))
    Phi <- I[LambdaSS, ] %*% Phi %*% I[, LambdaSS]
    
  } # END if ( !is.null(Phi) ) 
  
  ## Return a list of the output
  list(Lambda = Lambda,
       PhiMat = Phi)
  
}  #END orderFactors
