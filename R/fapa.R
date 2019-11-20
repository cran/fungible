#' Iterated Principal Axis Factor Analysis (fapa)
#'
#' This function applies the iterated principal axis factoring method to extract an unrotated factor structure matrix.
#'
#' @param R (Matrix) A correlation matrix to be analyzed.
#' @param numFactors (Numeric) The number of factors to extract.
#' @param epsilon (Numeric) A numeric threshold to designate whether the function has converged. The default value is 1e-4.
#' @param communality (Character) The routine requires an initial estimate of the communality values. There are three options (see below) with "SMC" (i.e., squared multiple correlation) being the default.
#' \itemize{
#'   \item \strong{"SMC"}: Initial communalities are estimated by taking the squared multiple correlations of each indicator after regressing the indicator on the remaining variables. The following equation is employed to find the squared multiple correlation: \eqn{1 - 1 / diag(R^-1)}.
#'   \item \strong{"maxr"}: Initial communalities equal the largest (absolute value) correlation in each column of the correlation matrix.
#'   \item \strong{"unity"}: Initial communalities equal 1.0 for all variables.
#'}
#' @param maxItr (Numeric) The maximum number of iterations to reach convergence. The default is 15,000.
#'
#' @details
#' \itemize{
#'   \item \strong{Initial communality estimate}: The choice of the initial communality estimate can impact the resulting principal axis factor solution.
#'   \itemize{
#'     \item \strong{Impact on the Estimated Factor Structure}: According to Widaman and Herringer (1985), the initial communality estimate does not have much bearing on the resulting solution \emph{when a stringent convergence criterion is used}. In their analyses, a convergence criterion of .001 (i.e., slightly less stringent than the default of 1e-4) is sufficiently stringent to produce virtually identical communality estimates irrespective of the initial estimate used. Based on their findings, it is not recommended to use a convergence criterion lower than 1e-3.
#'     \item \strong{Impact on the Iteration Procedure}: The initial communality estimates have little impact on the \emph{final factor structure} but they can impact the iterated procedure. It is possible that poor communality estimates produce a non-positive definite correlation matrix (i.e., eigenvalues <= 0) whereas different communality estimates  result in a converged solution. If the fapa procedure fails to converge due to a non-positive definite matrix, try using different communality estimates before changing the convergence criterion.
#'   }
#'
#' }
#'
#' @return The main output is the matrix of unrotated factor loadings.
#' \itemize{
#'   \item \strong{loadings}: (Matrix) A matrix of unrotated factor loadings extracted via iterated principal axis factoring.
#'   \item \strong{h2}: (Vector) A vector containing the resulting communality values.
#'   \item \strong{iterations}: (Numeric) The number of iterations required to converge.
#'   \item \strong{converged}: (Logical) TRUE if the iterative procedure converged.
#'   \item \strong{faControl}: (List) A list of the control parameters used to generate the factor structure.
#'   \itemize{
#'     \item \strong{epsilon}: (Numeric) The convergence criterion used for evaluating each iteration.
#'     \item \strong{communality}: (Character) The method for estimating the initial communality values.
#'     \item \strong{maxItr}: (Numeric) The maximum number of allowed iterations to reach convergence.
#'   }
#' }
#'
#' @author
#' \itemize{
#'   \item Casey Giordano (Giord023@umn.edu)
#'   \item Niels G. Waller (nwaller@umn.edu)
#'}
#'
#' @family Factor Analysis Routines
#' 
#' @references Widaman, K. F., & Herringer, L. G. (1985). Iterative least squares estimates of communality: Initial estimate need not affect stabilized value. \emph{Psychometrika, 50}(4), 469-477.
#'
#' @examples
#' ## Generate an example factor structure matrix
#' lambda <- matrix(c(.62, .00, .00,
#'                    .54, .00, .00,
#'                    .41, .00, .00,
#'                    .00, .31, .00,
#'                    .00, .58, .00,
#'                    .00, .62, .00,
#'                    .00, .00, .38,
#'                    .00, .00, .43,
#'                    .00, .00, .37),
#'                  nrow = 9, ncol = 3, byrow = TRUE)
#'
#' ## Find the model implied correlation matrix
#' R <- lambda %*% t(lambda)
#' diag(R) <- 1
#'
#' ## Extract factors using the fapa function
#' Out1 <- fapa(R           = R,
#'              numFactors  = 3,
#'              communality = "SMC")
#'
#' ## Call fapa through the factExtract function
#' Out2 <- faX(R          = R,
#'             numFactors = 3,
#'             facMethod  = "fapa",
#'             faControl  = list(communality = "maxr",
#'                               epsilon     = 1e-4))
#'
#' ## Check for equivalence of the two results
#' all.equal(Out1$loadings, Out2$loadings)
#'
#'
#' @export

fapa <- function(R,
                 numFactors    = NULL,
                 epsilon       = 1e-4,
                 communality   = "SMC",
                 maxItr        = 15000) {
  
  ## ~~~~~~~~~~~~~~~~~~ ##
  #### Error Checking ####
  ## ~~~~~~~~~~~~~~~~~~ ##
  
  ## Specify number of factors to extract
  if ( is.null(numFactors) ) {
    stop("The 'numFactors' argument is not specified, please specify how many factors are to be extracted.")
  } # END if ( is.null(numFactors) )
  
  ## Check to ensure R is positive definite before doing solve
  eigenVal <- eigen(R)$value
  if ( communality == "SMC" && min(eigenVal) <= 1E-8 ) {
    warning("Inverting the correlation matrix for SMC communality estimates requires a positive-definite matrix.")
    communality <- "maxr"  # NGW March 4, 2019
  } # END if ( min(eigenVal) <= 1E-8 )
  
  ## Find initial communality estimate
  if (communality == "SMC")     hsq <- 1 - 1 / diag(solve(R))
  if (communality == "maxr")    hsq <- apply( abs((R - diag(ncol(R)))) , 2, max)
  if (communality == "unity")   hsq <- rep(1, ncol(R))
  
  ## ~~~~~~~~~~~~~~~~~~ ##
  #### Begin Function ####
  ## ~~~~~~~~~~~~~~~~~~ ##
  
  ## Initialize parameter for convergence evaluation
  error <- 1
  
  ## Initialize iteration parameter
  iter <- 1
  
  ## Iterated principal axis factoring
  while (error > epsilon && iter < maxItr) {
    
    
    ## Reduced correlation matrix with communalities on diagonal
    diag(R) <- hsq
    
    ## Find eigenvalues & eigenvectors for reduced R
    Eigs <- eigen(R)
    
    
    ## Cannot take diag() of a scalar, same analyses, slightly different format
    if (numFactors > 1) {
      
      
      # NGW March 12 2019
      # Set negative eigenvalues to zero
      EigVals <- Eigs$values
      EigVals[EigVals < 0] <- 0
      
      
      ## Lambda = A %*% D^1/2
      Lambda <- Eigs$vec[, 1:numFactors] %*% diag(sqrt(EigVals[1:numFactors]))
      
    } else {
      
      ## Linear algebra with scalars in R does not work.
      ## Same analysis, different notation
      Lambda <- matrix(Eigs$vec[, 1] * sqrt(Eigs$val[1]), ncol = 1)
      
    } # END if (numFactors > 1)
    
    ## Find communality estimate from the factor loadings
    NewHsq <- diag( tcrossprod(Lambda) )
    
    ## Change the way the error term is computed
    
    ## Determine difference between old and new hsq estimates
    error <- abs( sum(NewHsq) - sum(hsq) )
    
    ## Rewrite hsq for the next iteration
    hsq <- NewHsq
    
    ## Track number of iterations used
    iter <- iter + 1
    
  } # END while (error > epsilon && iter < maxItr) {
  
  
  ## NGW March 12 2019
  convergence <- FALSE
  if (error < epsilon && iter < maxItr) {
    ## Did the algorithm converge or hit its max iteration attempts?
    convergence <- TRUE
  }
  
  
  
  ## List of the final output
  list(loadings   = Lambda,
       h2         = hsq,
       iterations = iter,
       converged  = convergence,
       faControl  = list(epsilon     = epsilon,
                         communality = communality,
                         maxItr      = maxItr))
  
} # END fapa <- function(R, numFactors, epsilon, maxItr) {



