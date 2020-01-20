#' Inter-Battery Factor Analysis by the Method of Maximum Likelihood
#' 
#' This function conducts maximum likelihood inter-battery factor analysis using procedures described by Browne (1979). 
#' The unrotated solution can be rotated  (using the \pkg{GPArotation} package) 
#' from a user-specified number of random (orthogonal) starting configurations. 
#' Based on the resulting complexity function value, the function determines the 
#' number of local minima and, among these local solutions, will find the 
#' "global minimum" (i.e., the minimized complexity value from the finite 
#' number of solutions). See Details below for an elaboration on the global 
#' minimum. This function can also return bootstrap standard errors of the factor solution.
#'
#' @param X (Matrix) A raw data matrix (or data frame).
#' @param R (Matrix) A correlation matrix.
#' @param n (Numeric) Sample size associated with the correlation matrix. Defaults to n = NULL.
#' @param NVarX (Integer) Given batteries X and Y, \code{NVarX} denotes the number of variables in battery X. 
#' @param numFactors (Numeric) The number of factors to extract for subsequent rotation.
#' @param itemSort (Logical) if \code{itemSort = TRUE} the factor loadings will be sorted within batteries.
#' @param rotate (Character) Designate which rotation algorithm to apply. The 
#' following are available rotation options: "oblimin", "quartimin", "targetT", 
#' "targetQ", "oblimax", "entropy", "quartimax", "varimax", "simplimax", 
#' "bentlerT", "bentlerQ", "tandemI", "tandemII", "geominT", "geominQ", "cfT", 
#' "cfQ", "infomaxT", "infomaxQ", "mccammon", "bifactorT", "bifactorQ", and 
#' "none". Defaults to rotate = "oblimin". See \pkg{GPArotation} package for more 
#' details. Note that rotations ending in "T" and "Q" represent orthogonal and 
#' oblique rotations, respectively.
#' @param rotateControl (List) A list of control values to pass to the factor rotation algorithms.
#' \itemize{
#'   \item \strong{numberStarts}: (Numeric) The number of random (orthogonal) 
#'   starting configurations for the chosen rotation method (e.g., oblimin). The first
#'   rotation will always commence from the unrotated factors orientation.
#'   Defaults to numberStarts = 10. 
#'   \item \strong{gamma}: (Numeric) This is a tuning parameter (between 0 
#'   and 1, inclusive) for an oblimin rotation.  See the \pkg{GPArotation} 
#'   library's oblimin documentation for more details. Defaults to gamma = 0 
#'   (i.e., a quartimin rotation).
#'   \item \strong{delta}: (Numeric) This is a tuning parameter for the geomin
#'    rotation. It adds a small number (default = .01) to the squared factor 
#'    loadings before computing the geometric means in the discrepancy function.
#'   \item \strong{kappa}: (Numeric) The main parameterization of the 
#'   Crawford-Ferguson (CF) rotations (i.e., "cfT" and "cfQ" for orthogonal and 
#'   oblique CF rotation, respectively). Defaults to kappa = 0. 
#'   \item \strong{k}: (Numeric) A specific parameter of the simplimax rotation. 
#'   Defaults to k = the number of observed variables.
#'   \item \strong{standardize}: (Character) The standardization routine used 
#'   on the unrotated factor structure. The three options are "none", "Kaiser", 
#'   and "CM". Defaults to standardize = "none". 
#'   \itemize{
#'     \item \strong{"none"}: No standardization is applied to the unrotated 
#'     factor structure. 
#'     \item \strong{"Kaiser"}: Use a factor structure matrix that has been 
#'     normed by Kaiser's method (i.e., normalize all rows to have a unit length).
#'     \item \strong{"CM"}: Use a factor structure matrix that has been normed
#'      by the Cureton-Mulaik method.
#'   }
#'   \item \strong{epsilon}: (Numeric) The rotational convergence criterion to 
#'   use. Defaults to epsilon = 1e-5.
#'   \item \strong{power}: (Numeric) Raise factor loadings the the n-th power 
#'   in the \code{\link{promaxQ}} rotation. Defaults to power = 4.
#'   \item \strong{maxItr}: (Numeric) The maximum number of iterations for the 
#'   rotation algorithm. Defaults to maxItr = 15000.
#' }
#' @param Seed (Integer) Starting seed for the random number generator.
#' @inheritParams faMain
#' 
#' @details
#' \itemize{
#'   \item \strong{Global Minimum}: This function uses several random starting 
#'   configurations for factor rotations in an attempt to find the global 
#'   minimum solution. However, this function is not guaranteed to find the 
#'   global minimum. Furthermore, the global minimum solution need not be 
#'   more psychologically interpretable than any of the local solutions (cf. 
#'   Rozeboom, 1992). As is recommended, our function returns all local 
#'   solutions so users can make their own judgements.
#'   \item \strong{Finding clusters of local minima}: We find local-solution sets by sorting the rounded  
#'   rotation complexity values (to the number of  digits specified in the \code{epsilon} 
#'   argument of the \code{rotateControl} list) into sets with equivalent values. For example, 
#'   by default \code{epsilon = 1e-5.} and thus \code{} will only evaluate the complexity 
#'   values to five significant digits. Any differences beyond that value will not effect the final sorting. 
#' }
#'
#' @return The \code{faIB} function will produce abundant output in addition 
#' to the rotated inter-battery factor pattern and factor correlation matrices. 
#' \itemize{
#'   \item \strong{loadings}: (Matrix) The rotated inter-battery factor solution with the 
#'   lowest evaluated discrepancy function. This solution has the lowest 
#'   discrepancy function \emph{of the examined random starting configurations}. 
#'   It is not guaranteed to find the "true" global minimum. Note that multiple
#'    (or even all) local solutions can have the same discrepancy functions.
#'   \item \strong{Phi}: (Matrix) The factor correlations of the rotated factor 
#'   solution with the lowest evaluated discrepancy function (see Details).
#'   \item \strong{fit}: (Vector) A vector containing the following fit statistics:
#'   \itemize{
#'      \item \strong{chiSq}: Chi-square goodness of fit value (see Browne, 1979, for details). Note that we apply Lawley's (1959) correction when computing the chi-square value.
#'      \item \strong{DF}: Degrees of freedom for the estimated model. 
#'      \item \strong{p-value}: P-value associated with the above chi-square statistic.
#'      \item \strong{MAD}: Mean absolute difference between the model-implied and the sample across-battery correlation matrices. A lower value indicates better fit. 
#'      \item \strong{AIC}: Akaike's Information Criterion where a lower value indicates better fit. 
#'      \item \strong{BIC}: Bayesian Information Criterion where a lower value indicates better fit. 
#'   }
#'   \item \strong{R}: (Matrix) Returns the (possibly sorted) correlation matrix, useful when raw data are supplied. 
#'   If \code{itemSort = TRUE} then the returned matrix is sorted to be consistent with the factor loading matrix.
#'   \item \strong{Rhat}: (Matrix) The (possibly sorted) reproduced correlation matrix.If \code{itemSort = TRUE} then the returned matrix is sorted to be consistent with the factor loading matrix.
#'   \item \strong{Resid}: (Matrix) A (possibly sorted) residual matrix (R - Rhat) for the between battery correlations. 
#'   \item \strong{facIndeterminacy}: (Vector) A vector (with length equal to the number of factors)
#'   containing Guttman's (1955) index of factor indeterminacy for each factor. 
#'   \item \strong{localSolutions}: (List) A list containing all local solutions 
#'   in ascending order of their factor loadings, rotation complexity values (i.e., the first solution 
#'   is the "global" minimum). Each solution returns the 
#'   \itemize{
#'      \item \strong{loadings}: (Matrix) the factor loadings, 
#'      \item \strong{Phi}: (Matrix) factor correlations, 
#'      \item \strong{RotationComplexityValue}: (Numeric) the complexity value of the rotation algorithm, 
#'      \item \strong{facIndeterminacy}: (Vector) A vector of factor indeterminacy indices for each common factor, and 
#'      \item \strong{RotationConverged}: (Logical) convergence status of the rotation algorithm. 
#'      }
#'   \item \strong{numLocalSets} (Numeric) How many sets of local solutions
#'    with the same discrepancy value were obtained. 
#'   \item \strong{localSolutionSets}: (List) A list containing the sets of 
#'   unique local minima solutions. There is one list element for every unique 
#'   local solution that includes (a) the factor loadings matrix, (b) the factor 
#'   correlation matrix (if estimated), and (c) the discrepancy value of the rotation algorithm. 
#'   \item \strong{rotate} (Character) The chosen rotation algorithm.
#'   \item \strong{rotateControl}: (List) A list of the control parameters 
#'   passed to the rotation algorithm.
#'   \item \strong{unSpunSolution}: (List) A list of output parameters (e.g., loadings, Phi, etc) from 
#'   the rotated solution that was obtained by rotating directly from the unrotated (i.e., unspun) common factor orientation. 
#'   \item \strong{Call}: (call) A copy of the function call.
#' }
#' 
#' @references Boruch, R. F., Larkin, J. D., Wolins, L., & MacKinney, A. C. (1970). 
#' Alternative methods of analysis: Multitrait-multimethod data. \emph{Educational 
#' and Psychological Measurement, 30}(4), 833–853. 
#' https://doi.org/10.1177/0013164470030004055
#' @references Browne, M. W.  (1979).  The maximum-likelihood solution in inter-battery factor analysis. 
#' \emph{British Journal of Mathematical and Statistical Psychology, 32}(1), 75-86.  
#' @references Browne, M. W.  (1980).  Factor analysis of multiple batteries by maximum likelihood.  
#' \emph{British Journal of Mathematical and Statistical Psychology, 33}(2), 184-199.  
#' @references Browne, M. W. (2001). An overview of analytic rotation in 
#' exploratory factor analysis. \emph{Multivariate Behavioral Research, 36}(1), 111-150.
#' @references Burnham, K. P. & Anderson, D. R.  (2004).  Multimodel inference: Understanding AIC and BIC in model selection.  
#' \emph{Sociological methods and research, 33}, 261-304.  
#' @references Cudeck, R. (1982). Methods for estimating between-battery factors,
#' \emph{Multivariate Behavioral Research, 17}(1), 47-68. 10.1207/s15327906mbr1701_3
#' @references Cureton, E. E., & Mulaik, S. A. (1975). The weighted varimax 
#' rotation and the promax rotation. \emph{Psychometrika, 40}(2), 183-195.
#' @references Guttman, L. (1955). The determinacy of factor score matrices with 
#' implications for five other basic problems of common factor theory. 
#' \emph{British Journal of Statistical Psychology, 8}(2), 65-81.
#' @references Tucker, L. R.  (1958).  An inter-battery method of factor analysis.  
#' \emph{Psychometrika, 23}(2), 111-136.  

#' @family Factor Analysis Routines
#'
#' @author
#' \itemize{
#'   \item Niels G. Waller (nwaller@umn.edu)
#'   \item Casey Giordano (Giord023@umn.edu)
#'}
#'
#' @examples
#' 
#' # Example 1:
#' # Example from: Browne, M. W.  (1979). 
#' #
#' # Data originally reported in:
#' # Thurstone, L. L. & Thurstone, T. G. (1941). Factorial studies 
#' # of intelligence. Psychometric Monograph (2), Chicago: Univ. 
#' # Chicago Press.
#' 
#' R.XY <- matrix(c(
#'  1.00, .554, .227, .189, .461, .506, .408, .280, .241,
#'  .554, 1.00, .296, .219, .479, .530, .425, .311, .311,
#'  .227, .296, 1.00, .769, .237, .243, .304, .718, .730,
#'  .189, .219, .769, 1.00, .212, .226, .291, .681, .661,
#'  .461, .479, .237, .212, 1.00, .520, .514, .313, .245,
#'  .506, .530, .243, .226, .520, 1.00, .473, .348, .290,
#'  .408, .425, .304, .291, .514, .473, 1.00, .374, .306,
#'  .280, .311, .718, .681, .313, .348, .374, 1.00, .672,
#'  .241, .311, .730, .661, .245, .290, .306, .672, 1.00), 9, 9)
#'
#'
#' dimnames(R.XY) <- list(c( paste0("X", 1:4),
#'                          paste0("Y", 1:5)),
#'                        c( paste0("X", 1:4),
#'                          paste0("Y", 1:5)))
#'                          
#'     out <- faIB(R = R.XY,  
#'                 n = 710,
#'                 NVarX = 4, 
#'                 numFactors = 2,
#'                 itemSort = FALSE,
#'                 rotate = "oblimin",
#'                 rotateControl = list(standardize  = "Kaiser",
#'                                      numberStarts = 10),
#'                 Seed = 1)
#'
#'  # Compare with Browne 1979 Table 2.
#'  print(round(out$loadings, 2))
#'  cat("\n\n")
#'  print(round(out$Phi,2))
#'  cat("\n\n MAD = ", round(out$fit["MAD"], 2),"\n\n")
#'  print( round(out$facIndeterminacy,2) )
#'  
#'  
#'  # Example 2:
#'  ## Correlation values taken from Boruch et al.(1970) Table 2 (p. 838)
#'  ## See also, Cudeck (1982) Table 1 (p. 59)
#'  corValues <- c(
#'    1.0,
#'    .11,  1.0,
#'    .61,  .47, 1.0,
#'    .42, -.02, .18,  1.0,
#'    .75,  .33, .58,  .44, 1.0, 
#'    .82,  .01, .52,  .33, .68,  1.0,
#'    .77,  .32, .64,  .37, .80,  .65, 1.0,
#'    .15, -.02, .04,  .08, .12,  .11, .13, 1.0,
#'    -.04,  .22, .26, -.06, .07, -.10, .07, .09,  1.0,
#'    .13,  .21, .23,  .05, .07,  .06, .12, .64,  .40, 1.0,
#'    .01,  .04, .01,  .16, .05,  .07, .05, .41, -.10, .29, 1.0,
#'    .27,  .13, .18,  .17, .27,  .27, .27, .68,  .18, .47, .33, 1.0,
#'    .24,  .02, .12,  .12, .16,  .23, .18, .82,  .08, .55, .35, .76, 1.0,
#'    .20,  .18, .16,  .17, .22,  .11, .29, .69,  .20, .54, .34, .68, .68, 1.0)
#'  
#'  ## Generate empty correlation matrix
#'  BoruchCorr <- matrix(0, nrow = 14, ncol = 14)
#'  
#'  ## Add upper-triangle correlations
#'  BoruchCorr[upper.tri(BoruchCorr, diag = TRUE)] <- corValues
#'  BoruchCorr <- BoruchCorr + t(BoruchCorr) - diag(14)
#'  
#'  ## Add variable names to the correlation matrix
#'  varNames <- c("Consideration", "Structure", "Sup.Satisfaction", 
#'  "Job.Satisfaction", "Gen.Effectiveness", "Hum.Relations", "Leadership")
#'  
#'  ## Distinguish between rater X and rater Y
#'  varNames <- paste0(c(rep("X.", 7), rep("Y.", 7)), varNames)
#'  
#'  ## Add row/col names to correlation matrix
#'  dimnames(BoruchCorr) <- list(varNames, varNames)
#'  
#'  ## Estimate a model with one, two, and three factors
#'  for (jFactors in 1:3) {
#'    tempOutput <- faIB(R          = BoruchCorr,
#'                       n          = 111,
#'                       NVarX      = 7,
#'                       numFactors = jFactors,
#'                       rotate     = "oblimin",
#'                       rotateControl = list(standardize  = "Kaiser",
#'                                            numberStarts = 100))
#'    
#'    cat("\nNumber of inter-battery factors:", jFactors,"\n")
#'    print( round(tempOutput$fit,2) )
#'  } # END for (jFactors in 1:3) 
#'  
#'  ## Compare output with Cudeck (1982) Table 2 (p. 60)
#'  BoruchOutput <- 
#'    faIB(R             = BoruchCorr,
#'         n             = 111,
#'         NVarX         = 7,
#'         numFactors    = 2,
#'         rotate        = "oblimin",
#'         rotateControl = list(standardize = "Kaiser"))
#'  
#'  ## Print the inter-battery factor loadings
#'  print(round(BoruchOutput$loadings, 3)) 
#'  print(round(BoruchOutput$Phi, 3)) 
#'  
#'  
#'  
#' @import GPArotation
#' @import stats
#' @export

faIB  <-  
  function(X             = NULL,
           R             = NULL,
           n             = NULL,
           NVarX         = 4, 
           numFactors    = 2, 
           itemSort      = FALSE,
           rotate        = "oblimin",
           bootstrapSE   = FALSE,
           numBoot       = 1000,      ## Number of bootstrapped samples drawn
           CILevel       = .95,       ## Bootstrap SE confidence level
           rotateControl = NULL,
           Seed          = 1){

  ## Return a copy of the call function. 
  Call <- match.call()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #                     ----ERROR CHECKING ----
  
  ## Must give the raw data to perform bootstraps
  if ( bootstrapSE == TRUE && is.null(X)) {
    stop("The raw data are required to compute the bootstrapped standard errors.")
  } # END if ( bootstrapSE == TRUE && is.null(X))
  
  
  ## ____Check 'rotate' ----
  ## Is a correct rotation is specified
  if ( !is.character(rotate) ) {
    stop("The 'rotate' argument must be a character string. See ?faMain for available options.")
  } # END if ( !is.character(rotate) ) 
  
  ## Character string of all available rotate options
  PlausibleRotations <- c("oblimin",   "quartimin", "targetT",
                          "targetQ",   "oblimax",   "entropy",
                          "quartimax", "varimax",   "simplimax",
                          "bentlerT",  "bentlerQ",  "tandemI",
                          "tandemII",  "geominT",   "geominQ",
                          "promaxQ",   "cfT",       "cfQ",
                          "infomaxT",  "infomaxQ",  "mccammon",
                          "bifactorT", "bifactorQ", "none")
  
  ## If rotate argument is incorrectly specified, return error
  if ( rotate %in% PlausibleRotations == FALSE ) {
    stop("The 'rotate' argument is incorrectly specified. Check for spelling errors (case sensitive).")
  } # END if (rotate %in% PlausibleRotations == FALSE) {
  
  #----____END ERROR CHECKING----
  
  
  ## Set the seed for reproducibility
  set.seed(Seed)
  
  # Calculate R matrix if not given
  if(!is.null(X)) R <- cor(X)
  nVar <- ncol(R)
  
  
  
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
  ####             DEFINE  FUNCTIONS              ####
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
  
  ## Function to compute model fit indices
  ComputeModelFit <- function(R, 
                              NVarX, 
                              numFactors, 
                              Lambda, 
                              Phi){
    
    ## Define the number of variables in the model
    nVar <- nrow(Lambda)
    
    ## Compute model-implied (reduced) correlation matrix
    Rhat <- Lambda %*% Phi %*% t(Lambda)
    
    ## Find the residual correlation matrix
    Resid <- R - Rhat
    
    ## Residual XY correlations
    ResidRXY <- Resid[1:NVarX, (NVarX+1):nVar]
    
    ## Mean absolute difference of the residual XY correlation matrix
    MAD <- mean( abs(ResidRXY) )
    
    ## Redefine the number of X variables
    p <- NVarX
    
    ## Within X-battery correlation matrix
    R.XX <- R[1:p, 1:p]
    ## Within Y-Battery correlation matrix
    R.YY <- R[(p+1):nVar, (p+1):nVar]
    ## Between battery correlation matrix
    R.XY <- R[1:p, (p+1):nVar]
    
    
    # Squared canonical correlations
    rhoSq <- Re(eigen(R.XY %*% solve(R.YY) %*% t(R.XY) %*% solve(R.XX))$values)
    
    
    # Browne 1979 formula for corrected chi-square
    chiSqTmp <- 1
    for(i in (numFactors+1):NVarX ){
      chiSqTmp <- chiSqTmp * (1 - rhoSq[i])
    }
    
    # see Browne 1979 p. 84
    pstarInv <- 0
    for(i in 1:numFactors){
      pstarInv <- pstarInv + 1/rhoSq[i]
    }
    
    # Lawley 1959 correction for n
    ## NOTE: No correction if sample size (n) is not specified
    if ( !is.null(n) ) {
      
      nstar <- (n - 1 - numFactors - .5*(nVar+1) + pstarInv)
      chiSq <- -nstar * log(chiSqTmp)
      dfChiSq <- (NVarX - numFactors)*( (nVar - NVarX) - numFactors)
      
      if(dfChiSq <= 0) warning("\n\n*****Zero Degrees of Freedom*****\n")
      
      
      # check if rotation is orthogonal
      # orthogonal is a logical (T/F)
      orthogonal <- all.equal(Phi, diag(numFactors), 
                              check.attributes = FALSE, 
                              check.names = FALSE) == TRUE 
      
      
      # ~~~~~~~~~~~~~~AIC and BIC~~~~~~~~~~~~~~~~~ 
      # qk = parsimony corrections factor
      if (orthogonal) {
        qk <- numEstParams <- (nVar * numFactors)
      } else {
        qk <- numEstParams <- (nVar * numFactors) + numFactors*(numFactors - 1)/2
      } # END if (orthogonal) 
      
      #BIC
      BIC <- chiSq + log(n) * qk
      #AIC
      AIC <- chiSq + 2 * qk
      
    } else {
      ## If sample size is NULL, set the following as NA
      BIC <- AIC <- chiSq <- dfChiSq <- NA
    } # END if ( !is.null(n) ) 
    
    
    
    ## Return
    list(Rhat  = Rhat,
         Resid = Resid,
         MAD   = MAD,
         chiSq = chiSq,
         DF    = dfChiSq,
         AIC   = AIC,
         BIC   = BIC)
  }# ED Compute Model Fit
  
  
  ## Extract names of the indicators to rename final output
  GetVarNames <- function(){
    varNames <- NULL 
    ## If variable names are provided, retain them
    if ( !is.null(X) )          varNames <- colnames(X)
    if ( !is.null(R) )          varNames <- colnames(R)
    if(is.null(varNames )) varNames <- paste0("V", 1:nVar)
    
    varNames
  }   
  
  # Main code to compute MLE multiple battery factor loadings
  MBFact <- function(R, NVarX){
    NVar <- ncol(R)
    p    <- NVarX
    R.XX <- R[1:p, 1:p]
    R.YY <- R[(p+1):NVar, (p+1):NVar]
    R.XY <- R[1:p, (p+1):NVar]
    
    F <- solve(R.XX) %*% R.XY %*% solve(R.YY) %*% t(R.XY)
    VDV.F <- eigen(F)
    V.F <- VDV.F$vectors
    D.F <- VDV.F$values
    G <- diag(sqrt(diag(t(V.F) %*% R.XX %*% V.F)))
    B.X <- V.F %*% solve(G)
    B.Y <- solve(R.YY)%*% t(R.XY) %*% B.X %*% diag(1/sqrt(D.F))
    
    
    # Unrotated loadings for X (sqrt of sigular values)
    LX <-(R.XX %*% B.X)[, 1:numFactors, drop = FALSE] %*% 
      diag(D.F^.25)[1:numFactors, 1:numFactors]
    
    # Unrotated loadings for Y
    LY <- R.YY %*% B.Y[, 1:numFactors, drop = FALSE] %*% 
      diag(D.F^.25)[1:numFactors, 1:numFactors]
    
    
    # &&&& Why do we have to ask for the Real part of the complex number
    urLoadings <- Re(rbind(LX, LY))
    urLoadings
  }#END MBFact
  

  SortItems <- function(Lambda1, varNames, NVarX){
            soutX <- faSort(Lambda1[1:NVarX, ],
                            phi = NULL,
                            BiFactor = FALSE,
                            salient = 0.25, 
                            reflect = FALSE)
  
            Lambda1[1:NVarX, ] <- soutX$loadings
  
  
            soutY <- faSort(Lambda1[(NVarX + 1):nVar, ],
                            phi = NULL,
                            BiFactor = FALSE,
                            salient = 0.25, 
                            reflect = FALSE)
  
            Lambda1[(NVarX + 1):nVar, ] <- soutY$loadings
            
            newOrder <- c(soutX$sortOrder,
                          (soutY$sortOrder + NVarX))
  
            rownames(Lambda1) <- varNames[newOrder] 
  
             list( Lambda1  = Lambda1, 
                   newOrder = newOrder)
  }# END SortItems

  #----cnRotate----
  ## Assign the default values for the rotateControl list
  cnRotate <- list(numberStarts = 10,
                   gamma        = 0,
                   delta        = .01,
                   kappa        = 0,
                   k            = nVar,
                   standardize  = "none",
                   epsilon      = 1e-5,
                   power        = 4,
                   maxItr       = 15000)
  
  ## Test that rotateControl arguments are correctly specified
  if ( !is.null(rotateControl) ) {
    
    ## If rotateControl is specified & standardize is misspecified:
    if ( !is.null(rotateControl$standardize) &&
         rotateControl$standardize %in% c("none", "Kaiser", "CM") == FALSE) {
      
      ## Produce an error
      stop("The 'standardize' argument must be: 'none', 'Kaiser', or 'CM'.")
    } # END
    
    ## Number of names in the default rotateControl list
    cnLength <- length( names(cnRotate) )
    
    ## Number of all unique names across user rotateControl and default rotateControl lists
    allLength <- length( unique( c( names(rotateControl), names(cnRotate) ) ) )
    
    ## If the lengths differ, it is because of spelling issue
    if (cnLength != allLength) {
      
      ## Find the incorrect arg/s
      incorrectArgs <- which( (names(rotateControl) %in% names(cnRotate) ) == FALSE)
      
      # stop("The following arguments are not valid inputs for the list of rotateControl arguments: ", paste0( paste(c(names(rotateControl)[incorrectArgs]), collapse = ", "), "." ) )
      stop(paste("The following arguments are not valid inputs for the list of rotateControl arguments:", paste(c(names(rotateControl)[incorrectArgs]), collapse = ", "), collapse = ": " ) )
      
    } # END if (cnLength != allLength)
    
  } # END if ( !is.null(rotateControl) )
  
  ## Change the default values based on user-specified rotateControl arguments
  cnRotate[names(rotateControl)] <- rotateControl
  
  
  
  
  ## Compute Guttman's factor determinacy indices
  GuttmanIndices <- function(Lambda, 
                             PhiMat,
                             SampCorr) {
    ## Purpose: Compute Guttman (1955) factor indeterminacy indices
    ##
    ## Args:    Lambda: (Matrix) Rotated factor loadings matrix
    ##          PhiMat: (Matrix) Factor correlation matrix
    
    ## Fator structure (works for either oblique or orthogonal model)
    facStruct <- Lambda %*% PhiMat
    
    Rinv  <-try(solve(SampCorr), silent = TRUE)
    ## If non-invertible, return NA values instead of returning an error
    if ( any( class(Rinv) %in% "try-error") ) {
      # warning("\n\nEncountered a singular R matrix when computing factor indeterminancy values\n")
      return( rep(NA, ncol(Lambda) ) )
    } # END if ( any( class(Rinv) %in% "try-error") ) 
    
    ## Factor indeterminacy solution
    FI <- sqrt( diag( t(facStruct) %*% Rinv %*% facStruct))
    
    if(max(FI)>1) FI <- rep(NA, length(FI))
    FI
  } # END GuttmanIndices
  
  
  ## Generate a random orthonormal starting matrix
  randStart <- function(dimension) {
    qr.Q(qr(matrix(rnorm(dimension^2), dimension, dimension )))
  } # END randStart <- function(dimension) 
  
  
  ## Do the rotations, will be called multiple times
  internalRotate <- function(lambda, 
                             rotation, 
                             spinMatrix,
                             rotateControl) {
    ## Purpose: Simple rotation wrapper
    ##
    ## Args: lambda:        (matrix) unrotated loadings to rotate
    ##       rotation:      (character) which rotation to use
    ##       spinMatrix:    (Matrix) randomly spin unrotated factor loadings
    ##       rotateControl: (list) tuning parameters to pass to rotation
    ##
    
    ## Determine column dimensions
    matrixDim <- ncol(lambda)
    
    # ## Perform rotation with the specified parameters
    # ## if numberStarts = 1 then rotate from the unrotated
    # ## solution
    # if (rotateControl$numberStarts == 1) {
    #   RandomSpinMatrix <- diag(matrixDim)
    # } # END if (rotateControl$numberStarts == 1) 
    # 
    # if (rotateControl$numberStarts > 1) {
    #   RandomSpinMatrix <- randStart(matrixDim)
    # } # END if (rotateControl$numberStarts > 1) 
    
    switch(rotation,
           "none" = {
             list(loadings = lambda,
                  Phi      = diag( matrixDim ) ) ## Identity matrix
           },
           "oblimin" = {
             GPArotation::oblimin(lambda,
                                  Tmat      = spinMatrix,
                                  gam       = cnRotate$gamma,
                                  normalize = FALSE,
                                  eps       = cnRotate$epsilon,
                                  maxit     = cnRotate$maxItr)
           },
           "quartimin" = {
             GPArotation::quartimin(lambda,
                                    Tmat      = spinMatrix,
                                    maxit     = cnRotate$maxItr,
                                    eps       = cnRotate$epsilon,
                                    normalize = FALSE)
           },
           # "targetT" = {
           #   GPArotation::targetT(lambda,
           #                        Tmat      = spinMatrix,
           #                        Target    = targetMatrix,
           #                        normalize = FALSE,
           #                        eps       = cnRotate$epsilon,
           #                        maxit     = cnRotate$maxItr)
           # },
           # "targetQ" = {
           #   GPArotation::targetQ(lambda,
           #                        Tmat      = spinMatrix,
           #                        Target    = targetMatrix,
           #                        normalize = FALSE,
           #                        eps       = cnRotate$epsilon,
           #                        maxit     = cnRotate$maxItr)
           # },
           "oblimax" = {
             GPArotation::oblimax(lambda,
                                  Tmat      = spinMatrix,
                                  normalize = FALSE,
                                  eps       = cnRotate$epsilon,
                                  maxit     = cnRotate$maxItr)
           },
           "entropy" = {
             GPArotation::entropy(lambda,
                                  Tmat      = spinMatrix,
                                  normalize = FALSE,
                                  eps       = cnRotate$epsilon,
                                  maxit     = cnRotate$maxItr)
           },
           "quartimax" = {
             GPArotation::quartimax(lambda,
                                    Tmat      = spinMatrix,
                                    normalize = FALSE,
                                    eps       = cnRotate$epsilon,
                                    maxit     = cnRotate$maxItr)
           },
           "varimax" = {
             GPArotation::Varimax(lambda,
                                  Tmat      = spinMatrix,
                                  normalize = FALSE,
                                  eps       = cnRotate$epsilon,
                                  maxit     = cnRotate$maxItr)
           },
           "simplimax" = {
             GPArotation::simplimax(lambda,
                                    Tmat      = spinMatrix,
                                    normalize = FALSE,
                                    eps       = cnRotate$epsilon,
                                    k         = cnRotate$k,
                                    maxit     = cnRotate$maxItr)
           },
           "bentlerT" = {
             GPArotation::bentlerT(lambda,
                                   Tmat      = spinMatrix,
                                   maxit     = cnRotate$maxItr,
                                   eps       = cnRotate$epsilon,
                                   normalize = FALSE)
           },
           "bentlerQ" = {
             GPArotation::bentlerQ(lambda,
                                   Tmat      = spinMatrix,
                                   maxit     = cnRotate$maxItr,
                                   eps       = cnRotate$epsilon,
                                   normalize = FALSE)
           },
           "tandemI" = {
             GPArotation::tandemI(lambda,
                                  Tmat      = spinMatrix,
                                  maxit     = cnRotate$maxItr,
                                  eps       = cnRotate$epsilon,
                                  normalize = FALSE)
           },
           "tandemII" = {
             GPArotation::tandemII(lambda,
                                   Tmat      = spinMatrix,
                                   maxit     = cnRotate$maxItr,
                                   eps       = cnRotate$epsilon,
                                   normalize = FALSE)
           },
           "geominT" = {
             GPArotation::geominT(lambda,
                                  Tmat      = spinMatrix,
                                  delta     = cnRotate$delta,
                                  normalize = FALSE,
                                  eps       = cnRotate$epsilon,
                                  maxit     = cnRotate$maxItr)
           },
           "geominQ" = {
             GPArotation::geominQ(lambda,
                                  Tmat      = spinMatrix,
                                  delta     = cnRotate$delta,
                                  normalize = FALSE,
                                  eps       = cnRotate$epsilon,
                                  maxit     = cnRotate$maxItr)
           },
           "promaxQ" = {
             promaxQ(urLoadings  = lambda,
                     power       = cnRotate$power,
                     standardize = cnRotate$standardize,
                     epsilon     = cnRotate$epsilon,
                     maxItr      = cnRotate$maxItr)
           },
           "cfT" = {
             GPArotation::cfT(lambda,
                              Tmat      = spinMatrix,
                              kappa     = cnRotate$kappa,
                              maxit     = cnRotate$maxItr,
                              eps       = cnRotate$epsilon,
                              normalize = FALSE)
           },
           "cfQ" = {
             GPArotation::cfQ(lambda,
                              Tmat      = spinMatrix,
                              kappa     = cnRotate$kappa,
                              eps       = cnRotate$epsilon,
                              normalize = FALSE,
                              maxit     = cnRotate$maxItr)
           },
           "infomaxT" = {
             GPArotation::infomaxT(lambda,
                                   Tmat      = spinMatrix,
                                   normalize = FALSE,
                                   eps       = cnRotate$epsilon,
                                   maxit     = cnRotate$maxItr)
           },
           "infomaxQ" = {
             GPArotation::infomaxQ(lambda,
                                   Tmat      = spinMatrix,
                                   normalize = FALSE,
                                   eps       = cnRotate$epsilon,
                                   maxit     = cnRotate$maxItr)
           },
           "mccammon" = {
             GPArotation::mccammon(lambda,
                                   Tmat      = spinMatrix,
                                   normalize = FALSE,
                                   eps       = cnRotate$epsilon,
                                   maxit     = cnRotate$maxItr)
           },
           "bifactorT" = {
             GPArotation::bifactorT(lambda,
                                    Tmat      = spinMatrix,
                                    normalize = FALSE,
                                    eps       = cnRotate$epsilon,
                                    maxit     = cnRotate$maxItr)
           },
           "bifactorQ" = {
             GPArotation::bifactorQ(lambda,
                                    Tmat      = spinMatrix,
                                    normalize = FALSE,
                                    eps       = cnRotate$epsilon,
                                    maxit     = cnRotate$maxItr)
           })
  } # END internalRotate
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #             ----   MAIN PROGRAM----
  
  # ----____Get varNames----
  varNames <- GetVarNames()
  
  #----____Compute MB factor loadings----
  
  urLoadings <- MBFact(R, NVarX)
  
  
  
  ## If only 1 factor, rotate must equal "none"
  if ( numFactors == 1 ) rotate <- "none"
  
  ## ----____ Standardize loadings ####
  Stnd <- faStandardize(method = cnRotate$standardize,
                        lambda = urLoadings)
  
  ## Extract DvInv for later unstandardization
  lambda <- Stnd$lambda
  
  
  ##----____ Rotate loadings from random spins ####
  
  ## If only 1 factor, rotate must equal "none"
  if ( numFactors == 1 ) rotate <- "none"
  
  ## Pre-allocate a list for the different attempts
  starts <- vector("list",
                   cnRotate$numberStarts)
  
  ## Set the seed for reproducibility in random spin matrices
  set.seed(seed = Seed)
  
  ## Create a list of matrices to randomly spin the factor structure
  starts <- lapply(starts, function(x) randStart(dimension = ncol(lambda)))
  
  ## First start is always from the unrotated factor orientation
  starts[[1]] <- diag(ncol(lambda))
  
  
  ## For each list element, do the specified rotate and save all output
  starts <- lapply(starts, function(randSpinMatrix) {
    
    rotatedLambda <- internalRotate(lambda    = lambda,
                                    rotation      = rotate,
                                    spinMatrix    = randSpinMatrix,
                                    rotateControl = cnRotate)
    
    ## If an orthogonal model, make Phi identity matrix
    if ( is.null(rotatedLambda$Phi) ) {
      
      ## Identity matrix for Phi
      rotatedLambda$Phi <- diag(numFactors)
      
    } # END if ( is.null(rotatedLambda$Phi) ) 
    
    ## Unstandardize the estimated loadings
    rotatedLambda$loadings[] <- Stnd$DvInv %*% rotatedLambda$loadings[]
    
    ## Return whole list, not just loadings
    rotatedLambda
    
  }) # END lapply(starts, function(x) )
  
  ## ----____Order Unique Solutions ####
  
  ## Save the minimized Complexity function for each attempt
  ##    $Table[,2] is the value of criterion at each iteration, which.min
  ##    finds the minimum criterion value across each attempt
  
  ## Evaluate the complexity functions to find smallest value
  if (rotate != "none") {
    
    ## For each random start, find the evaluated Complexity function
    ComplexityFunc <- sapply(starts, function(attempt) min(attempt$Table[, 2]))
    
  } # if (rotate != "none") 
  
  ## No Complexity function to evaluate when rotate = "none"
  if (rotate == "none") {
    
    ComplexityFunc <- rep(1, cnRotate$numberStarts)
    
  } # END if (rotate == "none") {
  
  ## Find the sort order (minimum first) of the rotation 
  ## complexity fnc values
  
  ## Order func creates an ordered vector with scalars representing 
  ## pre-ordered location
  
  ## E.g., if the smallest value is originally the 3rd element,
  ## order will create a vector with the 3rd element as "1"
  sortedComplexityOrder <- order(ComplexityFunc,
                                 decreasing = FALSE)
  
  
  ## UnSpun is *always* 1st rotation. Find where "1" is in ordered vector
  UnSpunPosition <- which.min(sortedComplexityOrder)
  
  ## Create new list to hold the rotated factor output 
  ## (loadings, phi, complexity fnc etc)
  uniqueSolutions <- vector("list", cnRotate$numberStarts)
  
  ## Create list of output from local solutions
  for ( iternumStart in 1:cnRotate$numberStarts ){
    
    ## Determine which element of sortedComplexityOrder to grab
    ## Start with lowest Complexity value, end with highest
    num <- sortedComplexityOrder[iternumStart]
    
    ## Extract the relevant factor loadings and Phi matrices
    ## "starts" is the unsorted list of rotated output
    selected <- starts[[num]]
    
    # December 9, 2019 Note that order factors (Factor Sort) is not yet implimented
    #### -------____ Sort rotated factor solutions ####
    
    ## Select the factor loadings
    uniqueSolutions[[iternumStart]]$loadings <- selected$loadings
    
    ## Select factor correlations
    uniqueSolutions[[iternumStart]]$Phi <- selected$Phi
    
    ## Select complexity values
    uniqueSolutions[[iternumStart]]$RotationComplexityValue <- ComplexityFunc[num]
    
    ## ----____Compute Factor Indeterminancy ####
    ## if the user provides urloadings then do not 
    ## compute factor indeterminancy values
    
    
    ## Guttman's factor indeterminacy indices
    uniqueSolutions[[iternumStart]]$facIndeterminacy <- 
      GuttmanIndices(Lambda = selected$loadings, 
                     PhiMat = selected$Phi,
                     SampCorr = R )   # perhaps change name of R matrix
    
    
    ## Did the local optima solution converge
    uniqueSolutions[[iternumStart]]$RotationConverged <- selected$convergence
    
  } # END for (iternumStart in 1:length(sortedDiscOrder))
  
  ## Extracted rotation complexity values from uniqueSolutions to create solution sets
  DisVal <- unlist(lapply(uniqueSolutions, function(x) x$RotationComplexityValue))
  
  ## Determine the number of digits for rotation convergence
  ## NOTE: subtract 2 because the leading zero and decimal count as characters
  ## NOTE: Convert scientific notation (e.g., 1e-5) to decimal form (works if non-scientific notation)
  numDigits <- nchar( format(cnRotate$epsilon, scientific = FALSE) ) - 2
  
  ## Round the discrepancy values to specified number of digits
  DisVal <- round(x      = DisVal,
                  digits = numDigits)
  
  ## Determine the sets of local solutions
  localMins <- unique(DisVal)
  
  # Classify solutions into UNIQUE groups
  nGrp <- length(localMins)
  localGrps <- vector("list", nGrp)
  names(localGrps) <- as.character(localMins)
  
  ## Sort the discrep functions into the local minima buckets
  for (local in 1:nGrp) {
    
    ## Determine which local solution belong in this set
    whichBelong <- which(DisVal == localMins[local])
    
    ## Select those matrices, put them in localGrps
    localGrps[[local]] <- uniqueSolutions[whichBelong]
    
  } # END for (local in 1:nGrp)
  
  ## Select the loading matrix with the lowest discrep value
  ## uniqueSolutions is ordered so 1 is minimum discrep value
  minLambda  <- uniqueSolutions[[1]]$loadings
  minPhi     <- uniqueSolutions[[1]]$Phi
  ## CG EDITS (30 sept 19): Changed "facIndeter" to a column vector
  facIndeter <- data.frame("FacIndeterminacy" = uniqueSolutions[[1]]$facIndeterminacy)
  
  ## If factor indeter. not computed, give NA values
  if ( is.null(facIndeter) ) facIndeter <- rep(NA, ncol(minLambda))
  
  
  
  Lambda1 <- minLambda
  
  ## ----____Reflect factors ----  
  
  # If reflect = TRUE then reflect factors s.t.
  # salient loadings are positive
  
  ## Determine whether factors are negatively oriented
  
  if (numFactors > 1) {
    Dsgn <- diag(sign(colSums(Lambda1^3))) 
  } else {
    Dsgn <- matrix(sign(colSums(Lambda1^3)), 1, 1) 
  } # END if (numFactors > 1)  
  
  ## If factors negatively oriented, multiply by -1, else multiply by 1
  Lambda1 <- Lambda1 %*% Dsgn
  
  if (!is.null(minPhi)) {
    ## If factor is negative, reverse corresponding factor correlations
    minPhi <- Dsgn %*% minPhi %*% Dsgn
  } # END if (!is.null(minPhi)) 
  
  
  
  
  
  
  
  # ----____Sort Items?----
  # if no item sort, newOrder is original order  
  newOrder <- 1:nVar
  if(itemSort == TRUE){
     if(numFactors == 1){
        newOrderX <- sort.list(Lambda1[1:NVarX])
        newOrderY <- sort.list(Lambda1[(NVarX+1):nVar])
        newOrder <- c(newOrderX , (newOrderY+NVarX) )
        Lambda1 <- matrix(Lambda1[newOrder], nVar, 1)
        rownames(Lambda1) <- varNames[newOrder]
        colnames(Lambda1) <- "F1"
        PhiNames <- colnames(Lambda1) <- paste0("F", 1:ncol(Lambda1))
        dimnames(minPhi) <- list(PhiNames, PhiNames)
     } 
     if(numFactors > 1){   
        sorted_out <- SortItems(Lambda1, varNames, NVarX)
        Lambda1 <- sorted_out$Lambda1
        newOrder <- sorted_out$newOrder
        PhiNames <- colnames(Lambda1) <- paste0("F", 1:ncol(Lambda1))
        dimnames(minPhi) <- list(PhiNames, PhiNames)
     }
    }else{
      rownames(Lambda1) <- varNames
      PhiNames <- colnames(Lambda1) <- paste0("F", 1:ncol(Lambda1))
      dimnames(minPhi) <- list(PhiNames, PhiNames)
    } #END Sort Items?

  
  
  # if itemSort = TRUE the reOrder items in R
    R <- R[newOrder, newOrder]


  #---- ____Compute Model Fit----   
    fitOut <- ComputeModelFit(R = R, 
                              NVarX = NVarX, 
                              numFactors = numFactors, 
                              Lambda = Lambda1, 
                              Phi = minPhi)
  
      Rhat  <- fitOut$Rhat
      Resid <- fitOut$Resid
      MAD   <- fitOut$MAD
      chiSq <- fitOut$chiSq
      DF    <- fitOut$DF
      AIC   <- fitOut$AIC
      BIC   <- fitOut$BIC

  
  
  #### -------- BOOTSTRAP SETUP -------- ####
  
  ## If true, compute bootstrap standard errors
  if (bootstrapSE == TRUE) {
    
    #This reorders the cols of X if itemSort = TRUE
    X <- X[, newOrder]
    # Number of subjects
    nSubj <- nrow(X)
    
    # Lists/variables to hold output
    ## Hold factor loadings from bootstraps
    fList <-
      ## Hold factor correlations from bootstraps
      phiList <- 
      ## Hold factor indeterminacy indices from bootstraps
      FIList <- vector("list", numBoot)
    
    ## Number of rows, used below for bootstrap sampling
    rows <- 1:nSubj
    
    #### ----____ Bootstrap for-loop ####
    
    ## Analyses on 'numBoot' number of random resamples of X matrix
    for (iSample in seq_len(numBoot)) {
      
      ## Set the seed for reproducibility
      set.seed(iSample + Seed)
      
      ## Resample (with replacement) from X (raw data matrix)
      bsSample <- sample(x       = rows,
                         size    = nSubj,
                         replace = TRUE)
      
      ## Create correlation matrix from the bootstrap sample
      Rsamp <- cor(X[bsSample, ])
      
      ## ____ Extract factor loadings ####
      
      ## Extract unrotated factors using resampled data matrix
      bsLambda <- MBFact(R     = Rsamp,
                         NVarX = NVarX)
      
      ## Conduct standardization
      bsStnd <- faStandardize(method = cnRotate$standardize,
                              lambda = bsLambda)
      
      ## Extract the standardized bootstrapped (unrotated) factor loadings
      bsLambda <- bsStnd$lambda
      
      ## ____ Rotate from random spins ####
      
      ## Find the "global" min out of all the random starting configurations
      
      ## Initialize a list to store rotations from random spins 
      bsStarts <- vector("list", cnRotate$numberStarts)
      
      ## Populate the list with random orthogonal matrices (i.e., start values for rotation)
      bsStarts <- lapply(bsStarts, function(x) randStart(dimension = numFactors))
      
      ## Conduct rotations from random start values
      bsStarts <- lapply(bsStarts, function(randSpinMatrix) {
        
        ## Rotate the bootstrapped samples
        bsRotated <- internalRotate(lambda        = bsLambda,
                                    rotation      = rotate,
                                    spinMatrix    = randSpinMatrix,
                                    rotateControl = cnRotate)
        
        ## If an orthogonal model, turn Phi into an identity matrix
        if ( is.null(bsRotated$Phi) ) {
          
          ## Identity matrix
          bsRotated$Phi <- diag(numFactors)
          
        } # END if ( is.null(bsRotated$Phi) ) 
        
        ## Return all output
        bsRotated
        
      }) # END bsStarts <- lapply(bsStarts, function(x)
      
      ## ____ Order unique solutions ####
      
      ## Find minimum of bsSolutions
      ## Evaluate the minimum disc functions to find smallest value
      if ( rotate != "none" ) {
        
        ## For each random start, find the evaluated discrepancy function
        bootstrapComplexityFunc <- sapply(bsStarts, function(attempt) min(attempt$Table[, 2]))
        
      }# END if ( rotate != "none" ) 
      
      ## If no rotation, no discrepancy value to evaluate
      if (rotate == "none") bootstrapComplexityFunc <- rep(1, cnRotate$numberStarts)
      
      ## Of all random configs, determinine which has the lowest criterion value
      bsMinimum <- which.min(bootstrapComplexityFunc)
      
      ## Extract global minimum factor loadings and Phi matrices from 
      bsLambda <- bsStarts[[bsMinimum]]$loadings
      bsPhi    <- bsStarts[[bsMinimum]]$Phi
      
      ## Unstandardize the rotated solution for this bootstrap sample
      bsLambda <- (bsStnd$DvInv %*% bsLambda)
      
      ## Pre-allocate list
      Aligned <- list()
      
      ## Factor alignment does not work on 1 factor (and doesn't make sense)
      if (numFactors > 1) {
        
        ## Align bootstrap sample with global minimum found earlier
        Aligned <- faAlign(F1   = Lambda1,
                           F2   = bsLambda,
                           Phi2 = bsPhi)
        
        ## Extract aligned elements
        AlignedLambda <- Aligned$F2
        AlignedPhi    <- Aligned$Phi2
        
      } # END if (numFactors > 1) 
      
      ## Cannot align 1-factor model, but can reflect the factor
      if (numFactors == 1) {
        if ( sum(bsLambda^3) < 0 ) bsLambda <- bsLambda * -1
        
        ## Define newly aligned elements
        AlignedLambda <- bsLambda
        AlignedPhi    <- bsPhi
      } # END if (numFactors == 1)
      
      ## Save loadings as the bootstrapped factor loadings
      fList[[iSample]]   <- AlignedLambda
      
      ## Save correlations as the bootstrapped factor correlations
      phiList[[iSample]] <- AlignedPhi
      
      ## Save factor indeterminacy indices
      FIList[[iSample]] <- GuttmanIndices(Lambda   = AlignedLambda,
                                          PhiMat   = AlignedPhi,
                                          SampCorr = Rsamp)
      
    } # END for (iSample in seq_len(numBoot))
    
    #### -----____ BootStrap stnd errors ####
    
    ## Convert list of matrices into array
    loadArray <- array(unlist(fList), c(nVar, numFactors, numBoot))
    phiArray  <- array(unlist(phiList), c(numFactors, numFactors, numBoot))
    FIArray   <- array(unlist(FIList), c(1, numFactors, numBoot))
   
    
    ## Bootstrap standard errors for factor loadings
    loadSE <- apply(loadArray, 1:2, sd)
    
    # Bootstrap standard errors for factor correlations
    phiSE <- apply(phiArray, 1:2, sd)
  
    
    ## ____ Confidence intervals ####
    
    ## Set alpha level
    alpha <- 1 - CILevel
    
    ## Confidence interval matrices for factor loadings
    loadCI.upper <- apply(loadArray, 1:2, quantile, 1 - (alpha / 2))
    loadCI.lower <- apply(loadArray, 1:2, quantile, (alpha / 2))
    
    # Confidence interval matrices for factor correlations
    phiCI.upper <- apply(phiArray, 1:2, quantile, 1 - (alpha / 2))
    phiCI.lower <- apply(phiArray, 1:2, quantile, (alpha / 2))
    
    ## ____ Bootstrap Fac Indeterminacy SE ####
    
    ## If none of the fac indeter values are NA, proceed. Else, set as NA
    if (!any(is.na(FIArray))) {
      
      ## Compute the standard error of the bootstrap factor indetermancy values
      FISE <- apply(FIArray, 1:2, sd)
      
      ## Find the upper/lower bounds of the confidence interval
      FICI.upper <- apply(FIArray, 2, quantile, 1 - (alpha / 2))
      FICI.lower <- apply(FIArray, 2, quantile, (alpha / 2))
      
      ## Create data frame of fac indeterminacy estimates
      facIndeter <- data.frame(facIndeter,
                               t(FISE),
                               FICI.lower,
                               FICI.upper)
      
      ## Add column names to data frame of FI estimates
      colnames(facIndeter) <- c("FI",
                                "SE",
                                paste0((alpha / 2) * 100, "th percentile"),
                                paste0((1 - (alpha / 2)) * 100, "th percentile")) 
      
      ## Add row names to correspond to each factor
      rownames(facIndeter) <-  paste0("f", 1:numFactors)
    }else{
      facIndeter <- NA
      FISE <- NA
    } # END if (!any(is.na(FIArray))) 
    
  
    #### ------- NAME BOOTSTRAP OBJECTS -------- ####
    
    ## Dimension names
    
    ## Name the factor indicators 
    rownames(loadSE) <-
      rownames(loadCI.upper) <-
      rownames(loadCI.lower) <- varNames[newOrder]
    
    ## Name the factors
    colnames(loadSE) <-
      colnames(loadCI.upper) <-
      colnames(loadCI.lower) <- 
      colnames(phiSE) <- rownames(phiSE) <-
      colnames(phiCI.upper) <- rownames(phiCI.upper) <-
      colnames(phiCI.lower) <- rownames(phiCI.lower) <-
      paste0("F", 1:numFactors)
    
  } # END if (bootstrapSE == TRUE)
  
  ## If not specified, return NULL
  if (bootstrapSE == FALSE) {
    loadSE          <-
      loadCI.upper  <-
      loadCI.lower  <-
      loadArray     <-
      phiSE         <-
      phiCI.upper   <-
      phiCI.lower   <-
      phiArray      <- 
      FIArray       <- 
      FISE          <- 
      FICI.upper    <- 
      FICI.lower    <- NULL
  } # END if (bootstrapSE == FALSE)
      
      
      ## CG Edits (13 dec 2019): If sample size is NULL, pchisq cannot be computed
      if ( !is.null(n) ) {
        fitVec <- c(chiSq, DF, pchisq(chiSq, DF), MAD, AIC, BIC)
      } else {
        fitVec <- c(chiSq, DF, NA, MAD, AIC, BIC)
      } # END if ( !is.null(n) ) 
      
      names(fitVec) <- c("chiSq", "DF", "p-value", "MAD", "AIC", "BIC")
      
      
      
  #----____Output ----   
  invisible(list(loadings              = Lambda1, # this is min complexity solutution
                 Phi                   = minPhi,
                 fit                   = fitVec,
                 R                     = R,
                 Rhat                  = Rhat,
                 Resid                 = Resid,
                 facIndeterminacy      = facIndeter,
                 loadingsSE            = loadSE,
                 CILevel               = CILevel,
                 loadingsCIupper       = loadCI.upper,
                 loadingsCIlower       = loadCI.lower,
                 PhiSE                 = phiSE,
                 PhiCIupper            = phiCI.upper,
                 PhiCIlower            = phiCI.lower,
                 facIndeterminacySE    = FISE,
                 localSolutions        = uniqueSolutions,
                 numLocalSets          = nGrp,
                 localSolutionSets     = localGrps,
                 rotate                = rotate,
                 loadingsArray         = loadArray,
                 PhiArray              = phiArray,
                 facIndeterminacyArray = FIArray,
                 rotateControl         = cnRotate,
                 unSpunSolution        = uniqueSolutions[[UnSpunPosition]],
                 Call = Call
     ))#END inivisble(list( . . .
  
} ##END faMB

