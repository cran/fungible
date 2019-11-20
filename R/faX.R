#' Factor Extraction (faX) Routines
#'
#' This function can be used to extract an unrotated factor structure matrix 
#' using the following algorithms: (a) unweighted least squares ("fals"); 
#' (b) maximum likelihood ("faml"); (c) iterated principal axis factoring ("fapa");
#' and (d) principal components analysis ("pca").
#'
#' @param R (Matrix) A correlation matrix used for factor extraction.
#' @param n (Numeric) Sample size associated with the correlation matrix. 
#' Defaults to n = NULL.
#' @param numFactors (Numeric) The number of factors to extract for subsequent rotation.
#' @param facMethod (Character) The method used for factor extraction. The 
#' supported options are "fals" for unweighted least squares, "faml" for maximum 
#' likelihood, "fapa" for iterated principal axis factoring, and "pca" for 
#' principal components analysis. The default method is "fals".
#' \itemize{
#'   \item \strong{"fals"}: Factors are extracted using the unweighted least 
#'   squares estimation procedure using the \code{\link{fals}} function.
#'   \item \strong{"faml"}: Factors are extracted using the maximum likelihood 
#'   estimation procedure using the \code{\link[stats]{factanal}} function.
#'   \item \strong{"faregLS"}: Factors are extracted using regularized 
#'   least squares factor analysis using the \code{\link{fareg}} function. 
#'   \item \strong{"faregML"}: Factors are extracted using regularized 
#'   maximum likelihood factor using the \code{\link{fareg}} function. 
#'   \item \strong{"fapa"}: Factors are extracted using the iterated principal 
#'   axis factoring estimation procedure using the \code{\link{fapa}} function.
#'   \item \strong{"pca"}: Principal components are extracted. 
#' }
#' 
#' @inheritParams faMain
#'
#'
#' @details
#' \itemize{
#'   \item \strong{Initial communality estimate}: According to Widaman and 
#'   Herringer (1985), the initial communality estimate does not have much 
#'   bearing on the resulting solution \emph{when the a stringent convergence 
#'   criterion is used}. In their analyses, a convergence criterion of .001 
#'   (i.e., slightly less stringent than the default of 1e-4) is sufficiently 
#'   stringent to produce virtually identical communality estimates irrespective 
#'   of the initial estimate used. It should be noted that all four methods for 
#'   estimating the initial communality in Widaman and Herringer (1985) are the 
#'   exact same used in this function. Based on their findings, it is not 
#'   recommended to use a convergence criterion lower than 1e-3.
#' }
#'
#' @return This function returns a list of output relating to the extracted factor loadings.
#' \itemize{
#'   \item \strong{loadings}: (Matrix) An unrotated factor structure matrix.
#'   \item \strong{h2}: (Vector) Vector of final communality estimates.
#'   \item \strong{faFit}: (List) A list of additional factor extraction output.
#'   \itemize{
#'     \item \strong{facMethod}: (Character) The factor extraction routine.
#'     \item \strong{df}: (Numeric) Degrees of Freedom from the maximum 
#'     likelihood factor extraction routine.
#'     \item \strong{n}: (Numeric) Sample size associated with the correlation matrix.
#'     \item \strong{objectiveFunc}: (Numeric) The evaluated objective function for the 
#'     maximum likelihood factor extraction routine. 
#'     \item \strong{RMSEA}: (Numeric) Root mean squared error of approximation 
#'     from Steiger & Lind (1980). Note that bias correction is computed if the 
#'     sample size is provided.
#'     \item \strong{testStat}: (Numeric) The significance test statistic for the maximum 
#'     likelihood procedure. Cannot be computed unless a sample size is provided. 
#'     \item \strong{pValue}: (Numeric) The p value associated with the significance test 
#'     statistic for the maximum likelihood procedure. Cannot be computed unless 
#'     a sample size is provided. 
#'     \item \strong{gradient}: (Matrix) The solution gradient for the least squares factor 
#'     extraction routine. 
#'     \item \strong{maxAbsGradient}: (Numeric) The maximum absolute value of the 
#'     gradient at the least squares solution. 
#'     \item \strong{Heywood}: (Logical) TRUE if a Heywood case was produced.
#'     \item \strong{converged}: (Logical) TRUE if the least squares or 
#'     principal axis factor extraction routine converged. 
#'   }
#' }
#'
#' @author
#' \itemize{
#'   \item Casey Giordano (Giord023@umn.edu)
#'   \item Niels G. Waller (nwaller@umn.edu)
#' }
#'
#' @family Factor Analysis Routines
#' 
#' @references Jung, S. & Takane, Y.  (2008).  Regularized common factor analysis.  
#' New trends in psychometrics, 141-149.  
#' @references Steiger, J. H., & Lind, J. (1980). Paper presented at the annual 
#' meeting of the Psychometric Society. \emph{Statistically-based tests for the 
#' number of common factors.}
#' @references Widaman, K. F., & Herringer, L. G. (1985). Iterative least squares 
#' estimates of communality: Initial estimate need not affect stabilized value. 
#' \emph{Psychometrika, 50}(4), 469-477.
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
#' ## Extract (principal axis) factors using the factExtract function
#' Out1 <- faX(R          = R,
#'             numFactors = 3,
#'             facMethod  = "fapa",
#'             faControl  = list(communality = "maxr",
#'                               epsilon     = 1e-4))
#'
#' ## Extract (least squares) factors using the factExtract function
#' Out2 <- faX(R          = R,
#'             numFactors = 3,
#'             facMethod  = "fals",
#'             faControl  = list(treatHeywood = TRUE))
#'
#' @export

faX <- function(R,
                n          = NULL,
                numFactors = NULL,
                facMethod  = "fals",
                faControl  = NULL){
               

  ## ~~~~~~~~~~~~~~~~~~ ##
  #### Error Checking ####
  ## ~~~~~~~~~~~~~~~~~~ ##
  
  ## Check R
  
  ## Symmetrical?
  if ( !isSymmetric(R) ) {
    stop("'R' is not a symmetric correlation matrix.")
  } # END if ( !isSymmetric(R) ) {
  
  
  ## Positive definite
  eigs <- eigen(R)$values
 
  
  if ( min(eigs) <= 0 && facMethod != "fals") {
    warning("R is not a positive definite matrix.")
  } # END if ( min(eigs) <= 0 )
  
  ## Check n
  
  if ( is.null(n) ) {
    n <- NA
  } # END if ( is.null(n) ) 
  
  ## Check facMethod
  
  ## Specified correctly?
  facMethodOptions <- c("fals", "faml", "fapa", "faregLS", "faregML", "pca")
  if ( (facMethod %in% facMethodOptions) == FALSE ) {
    stop("The method of factor extraction is incorrectly specified. Select either 'fals', 'faml', 'fapa', or 'pca'.")
  } # END if ( (facMethod %in% facMethodOptions) == FALSE )
  
  ## Check numFactors
  
  ## Specified?
  if ( is.null(numFactors) ) {
    stop("The 'numFactors' argument must be specified.")
  } # END if ( is.null(numFactors) )
  
  ## Check faControl
  
  ## Assert the default options
  cnFA <- list(treatHeywood   = TRUE,
               nStart         = 10,
               start          = NULL,
               maxCommunality = .995,
               epsilon        = 1e-4,
               communality    = "SMC",
               maxItr         = 15000)
  
  ## Correctly specified?
  if ( !is.null(faControl) ) {
    
    ## Total number of correct names
    cnFALength <- length( names(cnFA) )
    
    ## Total number of all names supplied
    allcnFALength <- length( unique( c( names(faControl), names(cnFA) ) ) )
    
    ## If lengths differ, something is mis-specified
    if (cnFALength != allcnFALength) {
      
      ## Find which are incorrect
      incorrectFAArgs <- which( (names(faControl) %in% names(cnFA)) == FALSE)
      
      stop(paste("The following arguments are not valid inputs for the list of factor analysis control arguments:", paste(c(names(faControl)[incorrectFAArgs]), collapse = ", "), collapse = ": " ) )
      
    } # END if (cnFALength != allcnFALength)
    
    ## Reassign values to the user-specified ones
    cnFA[names(faControl)] <- faControl
    
  } # END if ( !is.null(faControl) )
  
  

  
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
  #### Define Utility Func ####
  ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
  
  RMSEA <- function(discVal, DF, N = NULL) {
    ## Purpose: Compute root mean squared error of approximation
    ##          See Steiger & Lind (1980) & Browne & Cudeck (1993)
    ##
    ## Args:
    ##          discVal: discrepancy value from ML fac estimation
    ##          DF:      Degrees of freedom
    ##          N:       Scalar, sample size, used for bias correction
    ## 
    ## Output:
    ##          RMSEA: Scalar, RMSEA output
    
    ## If a sample size is listed, do bias correction
    if ( is.numeric(N) ) {
      ## Correct the discrepancy function for bias
      Discrepancy <- discVal - ( DF / N )
      
      ## Bias correct *can* produce negative values. Ensure value >= 0
      F0 <- max(Discrepancy, 0)
    } else {
      ## If no sample size, no bias-correction performed
      F0 <- discVal
    } # END if (is.numeric(N)) {
    
    RMSEA <- sqrt(F0 / DF)
  } # END SLRMSEA <- function(FAOutput, N = NULL) {
  
  ## ~~~~~~~~~~~~~~~~~~ ##
  #### Begin Function ####
  ## ~~~~~~~~~~~~~~~~~~ ##
  
  Out <- switch(facMethod,
                "fals" = {
                  fals(R            = R,
                       nfactors     = numFactors,
                       TreatHeywood = cnFA$treatHeywood)
                },
                "faml" = {
                  factanal(covmat   = R,
                           factors  = numFactors,
                           n.obs    = n,
                           start    = cnFA$start,
                           rotation = "none",
                           control  = list(nstart = cnFA$nStart,
                                           lower  = 1 - cnFA$maxCommunality))
                },
                 "faregLS" = {
                  fareg(R          = R,
                        numFactors = numFactors,
                        facMethod  = "rls")
                },
                "faregML" = {
                  fareg(R          = R,
                        numFactors = numFactors,
                        facMethod  = "rml")
                },
                "fapa" = {
                  fapa(R            = R,
                       numFactors   = numFactors,
                       epsilon      = cnFA$epsilon,
                       communality  = cnFA$communality,
                       maxItr       = cnFA$maxItr)
                }, 
                "pca" = {
                  
                  ## Extract eigen vectors and values
                  VLV <- eigen(R)
                  
                  ## Eivenvectors
                  # V <- VLV$vectors
                  V <- VLV$vectors[, seq_len(numFactors)]
                  
                  ## Eigenvalues
                  # L <- VLV$values
                  L <- VLV$values[seq_len(numFactors)]
                  
                  ## Find all principal components
                  loadings <- V %*% diag(sqrt(L))
                  
                  ## Compute communalities
                  h2 <- apply(loadings^2, 1, sum)
                  
                  ## Return output
                  list(loadings  = loadings,
                       h2        = h2,
                       converged = TRUE)
                })
  
  ## CG EDITS 27 Sept 2019: Re-update faControl after factor extraction
   ## Rationale: 'fapa' updates communality estimate if using SMC on NPD matrix
  cnFA[names(Out$faControl)] <- Out$faControl
  
  ## ------- Model Fit Indices -------- ## 
  
  modelFit <- list()
  
  ## Which factor extraction routine is used
  modelFit$facMethod <- facMethod
  
  ## Communalities
  modelFit$h2 <- apply(Out$loadings^2, 1, sum)
  
  ## Maximum likelihood indicies
  ## Degrees of freedom
  modelFit$df <- ifelse(facMethod == "faml", Out$dof, NA)
  ## Sample size
  modelFit$n  <- ifelse( !is.null(n), n, NA)
  ## Value of the objective function to be maximized
  modelFit$objectiveFunc <- ifelse(facMethod == "faml", 
                                   Out$criteria["objective"], NA)
  ## Steiger & Lind's root mean squared error of approximation
  modelFit$RMSEA <- ifelse(facMethod == "faml", 
                           RMSEA(discVal = modelFit$objectiveFunc,
                                 DF      = modelFit$df,
                                 N       = modelFit$n),
                           NA)
  ## ML test statistic
  modelFit$testStat <- ifelse(facMethod == "faml" && is.numeric(n), 
                              Out$STATISTIC, 
                              NA)
  ## ML test statistic's p value
  modelFit$pValue   <- ifelse(facMethod == "faml" && is.numeric(n), 
                              Out$PVAL, 
                              NA)
  
  ## Least squares indices
  
  ## Gradient 
  if (facMethod == "fals") {
    modelFit$gradient       <- Out$grdFALS
    modelFit$maxAbsGradient <- Out$MaxAbsGrad
  } else {
    modelFit$gradient       <- NA
    modelFit$maxAbsGradient <- NA
  } # END if (facMethod == "fals") 
  
  ## Were any Heywood cases detected
  modelFit$Heywood <- any(modelFit$h2 >= 1)
  
  ## Did the extraction procedure converge
  modelFit$converged <- Out$converged
  
  if (modelFit$converged == FALSE) {
    warning("The factor extraction method failed to converge.")
  } # END if (modelFit$converged == FALSE) 
  

  list(loadings = Out$loadings[],
       h2       = modelFit$h2,
       faFit    = modelFit,
       faControl = cnFA)
  
} # END factExtract

