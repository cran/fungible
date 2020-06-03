#' Multiple Battery Factor Analysis by Maximum Likelihood Methods
#' 
#' \code{faMB} estimates multiple battery factor analysis using maximum 
#' likelihood estimation procedures described by Browne (1979, 1980). Unrotated
#' multiple battery solutions are rotated (using the \pkg{GPArotation} package) 
#' from a user-specified number of of random (orthogonal) starting configurations. 
#' Based on procedures analogous to those in the \code{\link{faMain}} function,
#' rotation complexity values of all solutions are ordered to determine
#' the number of local solutions and the "global" minimum solution (i.e., the 
#' minimized rotation complexity value from the finite number of solutions). 
#'
#' @param X (Matrix) A raw data matrix (or data frame) structured in a subject 
#' (row) by variable (column) format. Defaults to \code{X = NULL}.
#' @param R (Matrix) A correlation matrix. Defaults to \code{R = NULL}.
#' @param n (Numeric) Sample size associated with either the raw data (X) or 
#' the correlation matrix (R). Defaults to \code{n = NULL}.
#' @param NB (Numeric) The number of batteries to analyze. In interbattery factor analysis NB = 2.
#' @param NVB (Vector) The number of variables in each battery. For example, 
#' analyzing three batteries including seven, four, and five variables 
#' (respectively) would be specified as \code{NVB = c(7, 4, 5)}.
#' @param numFactors (Numeric) The number of factors to extract for subsequent 
#' rotation. Defaults to \code{numFactors = NULL}.
#' @param epsilon (Numeric) The convergence threshold for the Gauss-Seidel iterator
#' when analyzing three or more batteries. Defaults to \code{epsilon = 1e-06}.
#' @param rotate (Character) Designate which rotation algorithm to apply. The 
#' following are available rotation options: "oblimin", "quartimin", 
#' "oblimax", "entropy", "quartimax", "varimax", "simplimax", 
#' "bentlerT", "bentlerQ", "tandemI", "tandemII", "geominT", "geominQ", "cfT", 
#' "cfQ", "infomaxT", "infomaxQ", "mccammon", "bifactorT", "bifactorQ", and 
#' "none". Defaults to rotate = "oblimin". See \pkg{GPArotation} package for more 
#' details. Note that rotations ending in "T" and "Q" represent orthogonal and 
#' oblique rotations, respectively.
#' @param PrintLevel (Numeric) When a value greater than zero is specified, 
#' \code{PrintLevel} prints the maximum change in communality estimates 
#' for each iteration of the Gauss-Seidel function. Note that Gauss-Seidel 
#' iteration is only called when three or more 
#' batteries are analyzed. Defaults to \code{PrintLevel = 0}.
#' @param Seed (Integer) Starting seed for the random number generator. 
#' Defaults to \code{Seed = 1}.
#'
#' @inheritParams faMain 
#'
#' @return The \code{faMB} function will produce abundant output in addition
#' to the rotated multiple battery factor pattern and factor correlation matrices. 
#' 
#' \itemize{
#'   \item \strong{loadings}: (Matrix) The (possibly) rotated multiple battery factor solution with the 
#'     lowest evaluated complexity value \emph{of the examined random starting configurations}. 
#'     It is not guaranteed to find the "true" global minimum. Note that multiple
#'      (or even all) local solutions can have the same discrepancy functions.
#'   \item \strong{Phi}: (Matrix) The factor correlations of the rotated factor 
#'     solution with the lowest evaluated discrepancy function (see Details).
#'   \item \strong{fit}: (Vector) A vector containing the following fit statistics:
#'   \itemize{
#'     \item \strong{ChiSq}: Chi-square goodness of fit value. 
#'     Note that, as recommended by Browne (1979), we apply Lawley's (1959) correction when computing the chi-square value when \code{NB = 2}.
#'     \item \strong{DF}: Degrees of freedom for the estimated model. 
#'     \item \strong{pvalue}: P-value associated with the above chi-square statistic.
#'     \item \strong{AIC}: Akaike's Information Criterion where a lower value indicates better fit. 
#'     \item \strong{BIC}: Bayesian Information Criterion where a lower value indicates better fit. 
#'     \item \strong{RMSEA}: Root mean squared error of approximation (Steiger & Lind, 1980).
#'   }
#'   \item \strong{R}: (Matrix) The \emph{sample} correlation matrix, 
#'   useful when raw data are supplied. 
#'   \item \strong{Rhat}: (Matrix) The \emph{reproduced} correlation matrix with communalities on the diagonal. 
#'   \item \strong{Resid}: (Matrix) A residual matrix (R - Rhat). 
#'   \item \strong{facIndeterminacy}: (Vector) A vector (with length equal to the number of factors)
#'   containing Guttman's (1955) index of factor indeterminacy for each factor. 
#'   \item \strong{localSolutions}: (List) A list (of length equal to the 
#'   \code{numberStarts} argument within \code{rotateControl}) containing all local solutions 
#'   in ascending order of their rotation complexity values (i.e., the first solution 
#'   is the "global" minimum). Each solution returns the following:
#'   \itemize{
#'      \item \strong{loadings}: (Matrix) the factor loadings, 
#'      \item \strong{Phi}: (Matrix) factor correlations, 
#'      \item \strong{RotationComplexityValue}: (Numeric) the complexity value of the rotation algorithm, 
#'      \item \strong{facIndeterminacy}: (Vector) A vector of factor indeterminacy indices for each common factor, and 
#'      \item \strong{RotationConverged}: (Logical) convergence status of the rotation algorithm. 
#'      }
#'   \item \strong{numLocalSets}: (Numeric) An integer indicating how many sets of local solutions
#'    with the same discrepancy value were obtained. 
#'   \item \strong{localSolutionSets}: (List) A list (of length equal to the 
#'   \code{numLocalSets}) that contains all local solutions with the same 
#'   rotation complexity value. Note that it is not guarenteed that all 
#'   solutions with the same complexity values have equivalent factor loading patterns. 
#'   \item \strong{rotate}: (Character) The chosen rotation algorithm.
#'   \item \strong{rotateControl}: (List) A list of the control parameters 
#'   passed to the rotation algorithm.
#'   \item \strong{unSpunSolution}: (List) A list of output parameters (e.g., loadings, Phi, etc) from 
#'   the rotated solution that was obtained by rotating directly from the unspun 
#'   (i.e., not multiplied by a random orthogonal transformation matrix) common 
#'   factor orientation. 
#'   \item \strong{Call}: (call) A copy of the function call.
#' 
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
#' @references Browne, M. and Cudeck, R. (1992). Alternative ways of assessing model fit. 
#' \emph{Sociological Methods and Research, 21(2)}, 230-258.
#' @references Burnham, K. P. & Anderson, D. R.  (2004).  Multimodel inference: Understanding AIC and BIC in model selection.  
#' \emph{Sociological methods and research, 33}, 261-304.  
#' @references Cudeck, R. (1982). Methods for estimating between-battery factors,
#' \emph{Multivariate Behavioral Research, 17}(1), 47-68. 10.1207/s15327906mbr1701_3
#' @references Cureton, E. E., & Mulaik, S. A. (1975). The weighted varimax 
#' rotation and the promax rotation. \emph{Psychometrika, 40}(2), 183-195.
#' @references Guttman, L. (1955). The determinacy of factor score matrices with 
#' implications for five other basic problems of common factor theory. 
#' \emph{British Journal of Statistical Psychology, 8}(2), 65-81.
#' @references Steiger, J. & Lind, J. (1980). Statistically based tests for the
#'  number of common factors. In \emph{Annual meeting of the Psychometric Society, 
#'  Iowa City, IA, volume 758}.
#' @references Tucker, L. R.  (1958).  An inter-battery method of factor analysis.  
#' \emph{Psychometrika, 23}(2), 111-136.  
#' 
#' @family Factor Analysis Routines
#' 
#' @author
#' \itemize{
#'   \item Niels G. Waller (nwaller@umn.edu)
#'   \item Casey Giordano (Giord023@umn.edu)
#'}
#' 
#' @examples
#' # These examples reproduce published multiple battery analyses. 
#' 
#' # ----EXAMPLE 1: Browne, M. W. (1979)----
#' #
#' # Data originally reported in:
#' # Thurstone, L. L. & Thurstone, T. G. (1941). Factorial studies 
#' # of intelligence. Psychometric Monograph (2), Chicago: Univ. 
#' # Chicago Press.
#' 
#' ## Load Thurstone & Thurstone's data used by Browne (1979)
#' data(Thurstone41)
#' 
#' Example1Output <-  faMB(R             = Thurstone41, 
#'                         n             = 710,
#'                         NB            = 2, 
#'                         NVB           = c(4,5), 
#'                         numFactors    = 2,
#'                         rotate        = "oblimin",
#'                         rotateControl = list(standardize = "Kaiser"))
#'                         
#' summary(Example1Output, PrintLevel = 2)                         
#' 
#' # ----EXAMPLE 2: Browne, M. W. (1980)----
#' # Data originally reported in:
#' # Jackson, D. N. & Singer, J. E. (1967). Judgments, items and 
#' # personality. Journal of Experimental Research in Personality, 20, 70-79.
#' 
#' ## Load Jackson and Singer's dataset
#' data(Jackson67)
#' 
#' 
#' 
#' Example2Output <-  faMB(R             = Jackson67, 
#'                         n             = 480,
#'                         NB            = 5, 
#'                         NVB           = rep(4,5), 
#'                         numFactors    = 4,
#'                         rotate        = "varimax",
#'                         rotateControl = list(standardize = "Kaiser"),
#'                         PrintLevel    = 1)
#'                         
#' summary(Example2Output)                         
#' 
#' 
#' 
#' # ----EXAMPLE 3: Cudeck (1982)----
#' # Data originally reported by:
#' # Malmi, R. A., Underwood, B. J., & Carroll, J. B. (1979).
#' # The interrelationships among some associative learning tasks. 
#' # Bulletin of the Psychonomic Society, 13(3), 121-123. DOI: 10.3758/BF03335032 
#' 
#' ## Load Malmi et al.'s dataset
#' data(Malmi79)
#' 
#' Example3Output <- faMB(R             = Malmi79, 
#'                        n             = 97,
#'                        NB            = 3, 
#'                        NVB           = c(3, 3, 6), 
#'                        numFactors    = 2,
#'                        rotate        = "oblimin",
#'                        rotateControl = list(standardize = "Kaiser"))
#'                        
#' summary(Example3Output)                        
#' 
#' 
#' 
#' # ----Example 4: Cudeck (1982)----
#' # Data originally reported by: 
#' # Boruch, R. F., Larkin, J. D., Wolins, L. and MacKinney, A. C. (1970). 
#' #  Alternative methods of analysis: Multitrait-multimethod data. Educational 
#' #  and Psychological Measurement, 30,833-853.
#' 
#' ## Load Boruch et al.'s dataset
#' data(Boruch70)
#' 
#' Example4Output <- faMB(R             = Boruch70,
#'                        n             = 111,
#'                        NB            = 2,
#'                        NVB           = c(7,7),
#'                        numFactors    = 2,
#'                        rotate        = "oblimin",
#'                        rotateControl = list(standardize  = "Kaiser",
#'                                             numberStarts = 100))
#'                                             
#' summary(Example4Output, digits = 3)                                             
#' 
#' @import GPArotation
#' @import stats
#' @export

faMB <- function(X = NULL,
                 R = NULL,
                 n = NULL,
                 NB = NULL, 
                 NVB = NULL, 
                 numFactors = NULL,
                 epsilon = 1E-06,
                 rotate = "oblimin",
                 rotateControl = NULL,
                 PrintLevel = 0,
                 Seed = 1){
  
  
  SampleSize = n
  NVar <- sum(NVB)
  NFac <- numFactors
  rotation <-rotate   

  
  ## Return a copy of the call function. 
  Call <- match.call()
  

  
 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #                     ----ERROR CHECKING ----
  
  ##----____Check X----
  
  ## Only check the following if X is specified
  if ( !is.null(X) ) {
    
    ## Missingness in data?
    if (nrow(X) != nrow(X[complete.cases(X),])) {
      warning("There are missing values in the data frame. See ?cor for how missing data is handled in the cor function.")
    } # END if (nrow(X) != nrow(X[complete.cases(X),]))
    
    ## Class must be matrix or DF to analyze via cor() function
    if ( !any( class(X) %in% c("matrix", "data.frame", "loadings") ) ) {
      stop("'X' must be of class matrix, data.frame, or loadings.")
    } # END if ( !any( class(X) %in% c("matrix", "data.frame", "loadings") ) )
    
  } # END if ( !is.null(X) )
  
  
  
  ## ----____Check R ----
  
  ## Make sure either corr mat or factor structure is supplied
  if ( is.null(R) && is.null(X) ) {
    stop("Must specify data to analyze using either the 'X'or  'R' arguments.")
  } # END if ( is.null(R) && is.null(X) )
  
  ## If R (corr matrix) is NOT supplied, make one using raw data
  if ( is.null(R) && !is.null(X) ) {
    
    ## Define correlation matrix
    R <- cor(X)
    
  } # END if ( is.null(R) && !is.null(X) )
  
  
  
  ## ----____Check variable names----
 
  ## Initialize variable in case no variable names are supplied
  VarNames <- NULL
  
  ## If X is supplied AND it has column names, use those
  if ( !is.null(X) && !is.null(colnames(X))) VarNames <- colnames(X)
  
  ## If R is supplied AND it has dimension names, use those
  if ( !is.null(R) && !is.null(dimnames(R))) {
    
    ## In case non-symmetric labels, pick dimension with label
    RowOrCol <- which.max( c(length(rownames(R)),
                             length(colnames(R))) )
    
    ## Define variable names based on whichever dimension(s) of R has labels
    VarNames <- dimnames(R)[[RowOrCol]]
    
    ## Set variable names in correlation matrix
    dimnames(R) <- list(VarNames, VarNames)
      
  } # END if ( !is.null(R) && !is.null(dimnames(R)))
  
  
  
  ## ----____Check NB ----
  
  ## IF number of batteries not specified, stop. 
  if ( is.null(NB) ) stop("\nMust specify the number of batteries in the 'NB' argument")
  
  
  
  ## ----____Check NVB ----
  
  ## IF number of batteries not specified, stop. 
  if ( is.null(NVB) ) stop("\nMust specify the number of variables per battery by the 'NVB' argument.")
  
  ## Ensure variables per battery argument is a vector
  if ( length(NVB) < 2) stop("\nMust correctly specify the number of variables per battery (i.e., NVB).")
  
  ## Length of NVB must equal NB
  if ( length(NVB) != NB) stop("\nThe number of variables per battery (NVB) must have the same number of elements as the number of batteries (NB).")
  
  
  ## ----____Check numFactors----
  
  ## Must specify numFactors
  if ( is.null(numFactors) ) {
    stop("The user needs to specify the number of factors to extract.")
  } # END if ( is.null(urLoadings) & is.null(numFactors) )
  
  if ( numFactors <= 0 ) {
    stop("\n\nFATAL ERROR:  Number of  factors must be greater than 0.")
  } # if ( numFactors <= 0 )
  
  
  # ----____Check 'rotate' ----
  
  ## If only 1 factor, rotate must equal "none"
  if ( numFactors == 1 ) rotate <- "none"
  
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
  
  
  


 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #           ---- DEFINE FUNCTIONS ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #   fncUpdatecnRotate       #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # ----UPDATE DEFAULT VALUES
  fncUpdatecnRotate <- function(rotateControl){
    ## Assign the default values for the rotateControl list
    cnRotate <- list(numberStarts = 10,
                     gamma        = 0,
                     delta        = .01,
                     kappa        = 0,
                     k            = NVar,
                     standardize  = "none",
                     epsilon      = 1e-5,
                     power        = 4,
                     maxItr       = 15000)
    
    
    ## Test that rotateControl arguments are correctly specified
    if ( !is.null(rotateControl) ){
      
      ## If rotateControl is specified & standardize is misspecified:
      if ( !is.null(rotateControl$standardize) &&
           rotateControl$standardize %in% c("none", "Kaiser", "CM") == FALSE){
        ## Produce an error
        stop("The 'standardize' argument must be: 'none', 'Kaiser', or 'CM'.")
      } # END if ( !is.null(rotateCo
      
      ## Number of names in the default rotateControl list
      cnLength <- length( names(cnRotate) )
      
      ## Number of all unique names across user rotateControl and default rotateControl lists
      allLength <- length( unique( c( names(rotateControl), names(cnRotate) ) ) )
      
      ## If the lengths differ, it is because of spelling issue
      if (cnLength != allLength){
        
        ## Find the incorrect arg/s
        incorrectArgs <- which( (names(rotateControl) %in% names(cnRotate) ) == FALSE)
        
        # stop("The following arguments are not valid inputs for the list of rotateControl arguments: ", paste0( paste(c(names(rotateControl)[incorrectArgs]), collapse = ", "), "." ) )
        stop(paste("The following arguments are not valid inputs for the list of rotateControl parameters:", paste(c(names(rotateControl)[incorrectArgs]), collapse = ", "), collapse = ": " ) )
        
      }# END if (cnLength != allLength)
      
    }#END if ( !is.null(rotateControl) )
    
    ## Change the default values based on user-specified rotateControl arguments
    
    cnRotate[names(rotateControl)] <- rotateControl
    
    #RETURN
    cnRotate
  }#END fncUpdatecnRotate
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #  fncRotateLocalSolutions  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  fncRotateLocalSolutions <- function(LambdaHat, cnRotate){
    ## If only 1 factor, rotate must equal "none"
    if ( numFactors == 1 ) rotate <- "none"
    
    ## ----Standardize loadings
    Stnd <- fungible::faStandardize(method = cnRotate$standardize,
                                    lambda = LambdaHat)
    
    ## Extract DvInv for later unstandardization
    LambdaHat <- Stnd$lambda
    
    ## Pre-allocate a list for the different attempts
    starts <- vector("list",
                     cnRotate$numberStarts)
    
    ## Set the seed for reproducibility in random spin matrices
    set.seed(seed = Seed)
    
    ## Create a list of matrices to randomly spin the factor structure
    starts <- lapply(starts, function(x) fncRandStart(dimension = numFactors))
    
    ## First start is always from the unrotated factor orientation
    starts[[1]] <- diag(numFactors)
    
    ## For each list element, do the specified rotate and save all output
    ListRotatedSolutions <- lapply(starts, function(randSpinMatrix){
      rotatedLambda <- fncRotate(lambda    = LambdaHat,
                                 rotation      = rotate,
                                 spinMatrix    = randSpinMatrix,
                                 rotateControl = cnRotate)
      
      ## If an orthogonal model, make Phi identity matrix
      if ( is.null(rotatedLambda$Phi) ){
        ## Identity matrix for Phi
        rotatedLambda$Phi <- diag(numFactors)
      } # END if ( is.null(rotatedLambda$Phi) ) 
      
      ## Unstandardize the estimated loadings
      rotatedLambda$loadings[] <- Stnd$DvInv %*% rotatedLambda$loadings[]
      
      ## Return whole list, not just loadings
      rotatedLambda
    }) # END lapply(starts, function(x) )
    
    #Return a list where each list element is a local solution 
    #(i.e., lambda, phi, facDetermin)
    ListRotatedSolutions
  }#END fncRotateLocalSolutions
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #  fncOrderUniqueSolutions  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  fncOrderUniqueSolutions <- function(ListRotatedSolutions, cnRotate){
    ## ----____Order Unique Solutions
    ## This function (1) orders the unique solutions
    ##               (2) Compute factor determinancy values for each solution
    ##               (3) Reflects factors for each solution 
    ##               (4) Create solutions sets with equal complexity values
    
    ## Save the minimized Complexity function for each attempt
    ##    $Table[,2] is the value of criterion at each iteration, which.min
    ##    finds the minimum criterion value across each attempt
    
    ## Evaluate the complexity functions to find smallest value
    if (rotate != "none") {
      
      ## For each random start, find the evaluated Complexity function
      ComplexityFunc <- sapply(ListRotatedSolutions, function(attempt) min(attempt$Table[, 2]))
      
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
    for( iternumStart in 1:cnRotate$numberStarts ){
      
      ## Determine which element of sortedComplexityOrder to grab
      ## Start with lowest Complexity value, end with highest
      num <- sortedComplexityOrder[iternumStart]
      
      ## Extract the relevant factor loadings and Phi matrices
      ## "ListRotatedSolutions" is the unsorted list of rotated output
      selected <- ListRotatedSolutions[[num]]
      
      # December 9, 2019 Note that order factors (Factor Sort) is not yet implimented
      ## Sort rotated factor solutions
      
      ## Select the factor loadings
      uniqueSolutions[[iternumStart]]$loadings <- selected$loadings
      
      ## Select factor correlations
      uniqueSolutions[[iternumStart]]$Phi <- selected$Phi
      
      ## Select complexity values
      uniqueSolutions[[iternumStart]]$RotationComplexityValue <- ComplexityFunc[num]
      
      ## Compute Factor Indeterminancy
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
    
  
    facIndeter <- data.frame("FacIndeterminacy" = uniqueSolutions[[1]]$facIndeterminacy)
    
    ## If factor indeter. not computed, give NA values
    if ( is.null(facIndeter) ) facIndeter <- rep(NA, ncol(minLambda))
    
    
    
    Lambda1 <- minLambda
    
    #Reflect Factors
    # Reflect factors s.t.
    # salient loadings are positive
    
    ## Determine whether factors are negatively oriented
    
    if (numFactors > 1){
      Dsgn <- diag(sign(colSums(Lambda1^3))) 
    }else{
      Dsgn <- matrix(sign(colSums(Lambda1^3)), 1, 1) 
    } # END if (numFactors > 1)  
    
    ## If factors negatively oriented, multiply by -1, else multiply by 1
    minLambda <- Lambda1 %*% Dsgn
    
    if (!is.null(minPhi)) {
      ## If factor is negative, reverse corresponding factor correlations
      minPhi <- Dsgn %*% minPhi %*% Dsgn
    } # END if (!is.null(minPhi)) 
    
    uniqueSolutions[[1]]$loadings <- minLambda
    uniqueSolutions[[1]]$Phi <- minPhi
    
    #fout <- uniqueSolutions[[1]]
    # Return
    list(uniqueSolutions = uniqueSolutions,
         nGrp = nGrp,
         localGrps = localGrps,
         UnSpunPosition = UnSpunPosition)
  }# END fncRotateLocalSolutions
  
  
  #~~~~~~~~~~~~~~~~~#
  # GuttmanIndices  #
  #~~~~~~~~~~~~~~~~~#
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
  
  
  
  #~~~~~~~~~~~~~~~~~#
  # fncRandStart    #
  #~~~~~~~~~~~~~~~~~#
  ## Generate a random orthonormal starting matrix
  fncRandStart <- function(dimension) {
    qr.Q(qr(matrix(rnorm(dimension^2), dimension, dimension )))
  } # END randStart <- function(dimension) 
  
  
  #~~~~~~~~~~~~~~~~~#
  # fncRotate       #
  #~~~~~~~~~~~~~~~~~#
  ## Do the rotations, will be called multiple times
   fncRotate  <- function(lambda, 
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
  } # END fncRotate
  
  
  #~~~~~~~~~~~~~~~~~#
  # fncBootstrap    #
  #~~~~~~~~~~~~~~~~~#
  fncBootstrap <- function(X, NB, NVar, NVB){
    
    S <- cor(X)
 
    # ----____Get Initial Estimates Lambda, V, mGamma, etc 
    initEst <- fncInitialEstimate(NB, NVB, NVar, S)
    Tmat <- initEst$Tmat
    mGamma <- initEst$mGamma
    V <- initEst$V
    
    
    # if NB = 2 this is MLE solution
    LambdaHat <- initEst$LambdaHat
    
    
    
    # ----Perform Gauss-Seidel Iterations
    if( NB > 2 ){
      
      # ----____Create Initial Mlist, Cij, Bij, 
      # needed for Chi square when NB = 2 and 
      # Gauss Seidel when NB > 2
      Mlist <- fncMlist(NB, NVB, mGamma, NFac) 
      
      # ----____Create BlockIJ
      BlockIJ <- fncBlockIJ(NB, NVB)  
      
      #Estimate LambdaHat 
      LambdaHatB <- fncGaussSeidel(Tmat, NVB, NB, V, Mlist, NFac, BlockIJ, PrintLevel)
      
      # Update Mlist a final time for computing X^2 values
      # mGamma is global from fncUpdateGamma function
      Mlist <- fncMlist(NB, NVB, mGamma, NFac) 
    } #END if (NB > 2)
    
    # return LambdaHat
    LambdaHatB
    
  }
  
  
  #~~~~~~~~~~~~~~~~~#
  # fncGaussSeidel  #
  #~~~~~~~~~~~~~~~~~#
  fncGaussSeidel <- function(Tmat, NVB, NB, V, Mlist, NFac, BlockIJ, epsilon, PrintLevel){
    
      # ----Estimate MLE loadings
      # ---- Gauss Seidel Iterations 
      NVar <- sum(NVB)
      comDelta <- comOld <-  rep(99, NVar)
      
     
       iter <- 1
       while( max(abs(comDelta)) > epsilon){
      
           # Update mGamma
           mGamma <- fncUpdateGamma(NFac, NB, V, Mlist, BlockIJ)
          
           #Update Mlist
           Mlist <- fncMlist(NB, NVB, mGamma, NFac)
          
           
           #TestCode
           TestCode <- F
           if(TestCode == TRUE){
              cat("\nmGamma iter", iter, "\n")
              print(mGamma[1:7, ])
              readline(prompt="\nPress [enter] to continue\n")
           
              cat("\nMlist iter", iter, "\n")
              print(Mlist[[2]])
              readline(prompt="\nPress [enter] to continue\n")
           }  
      
           # Assess convergence  
           communalities <- diag( mGamma %*% t(mGamma) )
           comDelta <- communalities - comOld
           comOld <- communalities
      
          if(PrintLevel > 0){
             cat(c(iter, "Max absolute change in communalities in iteraton", iter, "= ",  round(max(abs(comDelta)),6), "\n"))
           }
      
      iter <- iter + 1
      
     }# END while loop
    
    # ML estimate of unrotated loadings
    LambdaHat <- Tmat %*% mGamma 
    
    #return
    LambdaHat
    
  }#END   fncGaussSeidel
  

  
  #~~~~~~~~~~~~~~~~~#
  # fncChiSq2orLess #
  #~~~~~~~~~~~~~~~~~#
  fncChiSq2orLess <- function(R, 
                              NB,
                              NVB, 
                              numFactors, 
                              Lambda,
                              SampleSize){
    
      nVar <- sum(NVB)
      ## Within X-battery correlation matrix
      R.XX <- R[1:NVB[1], 1:NVB[1] ]
      ## Within Y-Battery correlation matrix
      R.YY <- R[(NVB[1]+1):nVar, (NVB[1]+1):nVar]
      ## Between battery correlation matrix
      R.XY <- R[1:NVB[1], (NVB[1]+1):nVar]
    
    
    # Squared canonical correlations
    rhoSq <- Re(eigen(R.XY %*% solve(R.YY) %*% t(R.XY) %*% solve(R.XX))$values)
    
    
    # Browne 1979 formula for corrected chi-square
    chiSqTmp <- 1
    for(i in (numFactors+1):NVB[1] ){
      chiSqTmp <- chiSqTmp * (1 - rhoSq[i])
    }
    
    # see Browne 1979 p. 84
    pstarInv <- 0
    for(i in 1:numFactors){
      pstarInv <- pstarInv + 1/rhoSq[i]
    }
    
    
    
    ChiSq <- dfChiSq <- pvalue <- NA
    
    # Lawley 1959 correction for n
    ## NOTE: No correction if sample size (n) is not specified
    if ( !is.null(SampleSize) ) {
      
      nstar <- (SampleSize - 1 - numFactors - .5*(nVar+1) + pstarInv)
      ChiSq <- -nstar * log(chiSqTmp)
      dfChiSq <- (NVB[1] - numFactors)*( (nVar - NVB[1]) - numFactors)
      pvalue <- pchisq(ChiSq, dfChiSq, lower.tail = FALSE)
      
      if(dfChiSq <= 0) warning("\n\n*****Zero Degrees of Freedom*****\n")
    }#END if ( !is.null(n) )
      
   
    
    ## Return
    list(ChiSq = ChiSq,
         DF    = dfChiSq,
         pvalue = pvalue, 
         rhoSq = rhoSq)
  }# END Compute Model Fit
  
  
  #~~~~~~~~~~~~~~~~~#
  # fncChiSq3orMore #
  #~~~~~~~~~~~~~~~~~#
  fncChiSq3orMore <- function(Mlist, NB, NVB, V, numFactors, SampleSize){
    
    NVar <- sum(NVB)
    Trace <- function(M) sum(diag(M))
    
    #structure of Mlist
    list("i", "Cij", "Bij", "Wij")
    
    #Browne 1980 EQ 3.26
    
    Cq <- lapply(Mlist, "[[", 2)
    
    # Create B, BVB
        Bq <- lapply(Mlist, "[[", 3)
        B <- matrix(0, 1, ncol(Bq[[1]]))
        iB <- length(Bq)

        for( iterB in 1:iB){
            B <- rbind(B, Bq[[iterB]])
        }
     B <- B[-1, ]
     BVB <- t(B) %*% V %*% B
    
    # Calculus Wsum
    Wq <- lapply(Mlist, "[[", 4)
    Imat <- diag(nrow(Wq[[1]]))
    

    DetIminusCq <- 1
    Wsum <- matrix(0, nrow = nrow(Wq[[1]]), ncol = ncol(Wq[[1]]))
   
     trWq <- 0
    for(i in 1: NB){ 
      Wsum <- Wsum + Wq[[i]]
      DetIminusCq <- DetIminusCq * det(Imat - Cq[[i]])
      trWq <- trWq + Trace(Wq[[i]])
    }
    
   F1 <- log( det(Imat + Wsum) * DetIminusCq  ) - log(det(V))
   F2 <- trWq - Trace(  BVB %*% solve(Imat + Wsum)  )
   
  
   
   Ffit <- F1 + F2
   
   ChiSq <- "Sample Size not provided"
   dfChiSq <- NA
   pvalue <- NA
   AIC <- NA
   BIC <- NA
   
   if(!is.null(SampleSize)){
       ChiSq <- (SampleSize - 1) * Ffit
       NVar <- sum(NVB)
       dfChiSq <- 1/2 * ( (NVar - numFactors)^2  - ( sum(NVB^2) + numFactors) )
       pvalue <- pchisq(ChiSq, dfChiSq, lower.tail = FALSE)
   }#END if(!is.null(SampleSize))
      
   
  ## Return
  list(ChiSq = ChiSq,
       DF    = dfChiSq,
       pvalue = pvalue,
       rhoSq = NA)
  
  }#END fncChiSq3orMore
  
  
  
  #~~~~~~~~~~~~~~~~~#
  # fncModelFit     #
  #~~~~~~~~~~~~~~~~~#
  # Function to compute model fit indices
  fncModelFit <- function(R, 
                          NB,
                          NVB,
                          Mlist, 
                          V,
                          numFactors, 
                          Lambda,
                          Phi,
                          SampleSize){
    
    if(NB <= 2){
         FitOut <- fncChiSq2orLess(R, 
                                   NB,
                                   NVB, 
                                   numFactors, 
                                   Lambda, 
                                   SampleSize)
    }
    if(NB >=3){

      FitOut <- fncChiSq3orMore(Mlist, 
                                NB,
                                NVB,
                                V,
                                numFactors, 
                                SampleSize)
    }  
    
      ChiSq <- FitOut$ChiSq
      dfChiSq    <- FitOut$DF
      pvalue <- FitOut$pvalue
      AIC   <- FitOut$AIC
      BIC   <- FitOut$BIC
      rhoSq <- FitOut$rhoSq
    
    ## Compute model-implied (reduced) correlation matrix
    # At this stage, Lambda is unrotated (i.e., orthogonal)    
    Rhat <- Lambda %*% Phi %*% t(Lambda)
    
    ## Find the residual correlation matrix
    Resid <- R - Rhat
  
    
    if ( !is.null(n) ) {
    
      # check if rotation is orthogonal
      # orthogonal is a logical (T/F)
      orthogonal <- all.equal(Phi, diag(numFactors), 
                              check.attributes = FALSE, 
                              check.names = FALSE) == TRUE 
      
      
      # ~~~~~~~~~~~~~~AIC and BIC~~~~~~~~~~~~~~~~~ 
      # qk = parsimony corrections factor
      if (orthogonal) {
        qk <- numEstParams <- (NVar * numFactors)
      } else {
        qk <- numEstParams <- (NVar * numFactors) + numFactors*(numFactors - 1)/2
      } # END if (orthogonal) 
      
      #BIC
      BIC <- ChiSq + log(n) * qk
      #AIC
      AIC <- ChiSq + 2 * qk
      
      
      ## RMSEA
      ## NOTE: If ChiSq < DF, set numerator to zero
      RMSEA <- sqrt( max(ChiSq - dfChiSq, 0) /
                         ((SampleSize -1 ) * dfChiSq ) )
      
    } else {
      ## If sample size is NULL, set the following as NA
      BIC <- AIC <- ChiSq <- dfChiSq <- RMSEA <- NA
    } # END if ( !is.null(n) ) 
    
    
    
    ## Return
    list(Rhat  = Rhat,
         Resid = Resid,
         ChiSq = ChiSq,
         DF    = dfChiSq,
         pvalue = pvalue,
         AIC   = AIC,
         BIC   = BIC,
         RMSEA = RMSEA,
         rhoSq = rhoSq)
  }# END fncModel Fit
  
  
  #~~~~~~~~~~~~~~~~~#
  # fncTmat         #
  #~~~~~~~~~~~~~~~~~#   
# fncTmat: Compute block diag of lower tri choleski matrices
  fncTmat <- function(NB, NVB, NVar, S){
    Tmat <- matrix(0, NVar, NVar)
    j <- 1
    k <- 0
    
    for(i in 1:NB){
        k <- k+NVB[i]
        Tmat[j:k, j:k] <- t(chol(S[j:k, j:k]))
        j <- j + NVB[i]
    } #end for loop
    
    
    Tmat
  }#END fncTmat

  #~~~~~~~~~~~~~~~~~#
  # fncBlockIJ      #
  #~~~~~~~~~~~~~~~~~# 
  fncBlockIJ <- function(NB, NVB){
    # Find start and end columns for batteries in supermatrix for 
    # computing V
    # See Cudeck 1982 EQ 12
    
    # matrix to hold supermatrix block indices
    BlockIJ <- matrix(0, NB^2, 6)
    
    iterK<- 1
    
    for(iterI in 1:NB){
      for(iterJ in 1: NB){
        
        BlockIJ[iterK, 1]    <- iterI
        BlockIJ[iterK, 2]    <- iterJ
        
        RowStart <- sum(NVB[1:iterI]) - NVB[iterI] + 1
        RowEnd    <- sum(NVB[1:iterI])
        
        ColStart <- sum(NVB[1:iterJ]) - NVB[iterJ] + 1
        ColEnd <-   sum(NVB[1:iterJ])
        
        BlockIJ[iterK, 3:4]  <- c(RowStart, RowEnd)
        BlockIJ[iterK, 5:6]  <- c(ColStart, ColEnd)
        
        iterK <- iterK + 1
      }#END for iterJ in 1:NB 
      
    } #END for iterI in 1:NB
    
    
    #rS = row Start index;    rE = row End index
    #cS = column Start index; cE = column End index
    
    colnames(BlockIJ) <- c("I", "J", 
                           "rS", "rE",
                           "cS", "cE")
    BlockIJ
  } #END  fncBlockIJ
 
   
  #~~~~~~~~~~~~~~~~~~~~#
  #     fncMlist       #
  #~~~~~~~~~~~~~~~~~~~~#
  fncMlist <- function(NB, NVB, mGamma, NFac){
    #Create list of LISTS OF MATRICES Cij, Bij, etc
    xlist <- list("i", "Cij", "Bij", "Wij")
    names(xlist) <- c("i", "Cij", "Bij", "Wij")
    
    Mlist <- list()
    
    Imat <- diag(NFac)
    
    for(iM in 1:NB){
       Mlist[[iM]] <- xlist
    } #END for(iM in 1:NB)
    
     j <- 1
     k <- 0
    
     for(iNB in 1:NB){
         k <- k + NVB[iNB]
         mGammaj <- mGamma[j:k, ]
         
         #check 
         #cat("\n j, k", j, k, "\n")

         Cj <- t(mGammaj) %*% mGammaj                   #Cudeck EQ 9
         InvImatMinusCj <- solve( Imat   - Cj )
         Bj <- mGammaj %*%  InvImatMinusCj              #Cudeck EQ 10
         Wj <- Cj %*% InvImatMinusCj                    #Cudeck EQ 11
      
         Mlist[[iNB]][[1]]  <- iNB
         Mlist[[iNB]][[2]]  <- Cj
         Mlist[[iNB]][[3]]  <- Bj
         Mlist[[iNB]][[4]]  <- Wj
      
         j <- j + NVB[iNB]
        }# END for(iNB in 1:NB)
    
    Mlist
  } # END fncMlist
  
  #~~~~~~~~~~~~~~~~~~~~#
  #     fncmakeWinv      #
  #~~~~~~~~~~~~~~~~~~~~#
  fncmakeWinv <- function(NFac, iBattery, Mlist){
    
      w <- lapply(Mlist, "[[", 4)
      Winv <- matrix(0, NFac, NFac)
      

      # make Winv
      for(jIterW in 1:NB){
          if(iBattery != jIterW) Winv = Winv + w[[jIterW]]
      }
  
    solve(Winv)
  }# END makeWinv 
  
  
  #~~~~~~~~~~~~~~~~~~~~#
  # fncUpdatemGamma    #
  #~~~~~~~~~~~~~~~~~~~~#
  fncUpdateGamma <- function(NFac, NB, V, Mlist, BlockIJ){  
    
    kRow <- 0
    
    for(iterBattery in 1:NB){
      
         #Initialize VBj
         VBj <- 0
      
      
         for(r in 1:NB){
             kRow <- kRow +1
        
             ## See Cudeck EQ 11 and 12 for update method look at 0's and 1's superscript
             if(r != iterBattery){
                 #Browne 1980 EQ 3.19
                VBj <- VBj + V[ (BlockIJ[kRow,3]: BlockIJ[kRow,4]), (BlockIJ[kRow,5]: BlockIJ[kRow,6])]  %*% Mlist[[r]]$Bij 
             } #END   if(r != iterBattery)
        
          }# END for(r in 1:NB)
      
      
      Winv <- fncmakeWinv(NFac, iBattery = iterBattery, Mlist)
      
      StartEnd <- which(BlockIJ[, 1] == iterBattery)[1] 
      
      
      ## January 6, 2020 For some reason updated mGamma needs to be global 
      ## (see double assignment below) for this
      ## code to work properly.  I do not know why
      mGamma[(BlockIJ[StartEnd,3]: BlockIJ[StartEnd,4]), ] <<- VBj %*% Winv
      
      Mlist <- fncMlist(NB, NVB, mGamma, NFac) 
     
      
    }#END for(iterBattery in 1:NB)
  
    #Return aGamma
    mGamma
  }# END  fncUpdateGamma
  
  
  #~~~~~~~~~~~~~~~~~~~~#
  # fncInitialEstimate #
  #~~~~~~~~~~~~~~~~~~~~#
  # Get initial estimates of V, mGamma and LambdaHat
  fncInitialEstimate <- function(NB, NVB, NVar, S){ 
    
  
    m <- NB
    
    # ---Generate T
    Tmat <- fncTmat(NB, NVB, NVar, S)
 
    
    TmatInv <- solve(Tmat)
    
    # ---- Create V
    V <- TmatInv %*% S %*% t(TmatInv)
    
    # Eigen decomposition of V
    UDU <- eigen(V)
    D <- diag(NFac)
    diag(D) <- UDU$values[1:NFac]
    I <- diag(nrow(D))
    U <- UDU$vectors
    
    
    # Create Initial mGamma
    mGamma <- as.numeric(sqrt(m/(m-1))) * U[,1:NFac] %*% sqrt(D-I)
    
    # Initial estimate of MB factor loadings
    # if NB ==2 this is the closed form MLE estimate
    LambdaHat <- Tmat %*% mGamma 
    
    
    list(Tmat = Tmat,
         V = V,
         mGamma = mGamma,
         LambdaHat = LambdaHat)
  } #END fncInitialEstimate
  

  #~~~~~~~~~~~~~~~~~~~~#
  # fncAddDimNames     #
  #~~~~~~~~~~~~~~~~~~~~#
  # Provide labels to the resulting output
  fncAddDimNames <- function(fout){
    
    ## If no user data has no variable names are supplied, remain unnamed
    VarNames <- NULL
    
    ## If R has row or column names, use those
      ## Note: prior error checking assures row/col have same labels
    if (!is.null(dimnames(R))) VarNames <- colnames(R)
    
    ## Default factor labels as "f1", "f2", etc.
    FacNames <- paste0("f", seq_len(numFactors))
    
    ## Apply variable and factor names to function output
    dimnames(fout$loadings) <- ## Factor loadings
      dimnames(fout$unSpunSolution$loadings) <- ## Unspun factor loadings
      list(VarNames, FacNames)
    dimnames(fout$Phi) <- ## Factor correlations
      dimnames(fout$unSpunSolution$Phi) <-  ## Unspun factor correlations
      list(FacNames, FacNames)
    dimnames(fout$Rhat) <- ## Model-implied correlation matrix
      dimnames(fout$Resid) <- ## Residual correlation matrix
      list(VarNames, VarNames)
    names(fout$facIndeterminacy) <- ## Factor indeterminacy
      FacNames
    
    ## Return
    fout
    
  } # END fncAddDimNames
  
  
  
  
  #~~~~~~~~~~ END FUNCTION DEFINITIONS~~~~~~~~~~~~~~
  
 
  
  
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
 #                    ----MAIN----
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
  
  # ----Update rotation control parameters----
  cnRotate <- fncUpdatecnRotate(rotateControl)
  
  # ----Get Initial Estimates Lambda, V, mGamma, etc ----
  initEst <- fncInitialEstimate(NB, NVB, NVar, S = R)
  Tmat <- initEst$Tmat
  mGamma <- initEst$mGamma
  V <- initEst$V
  # if NB = 2 this is the MLE solution
  LambdaHat <- initEst$LambdaHat
  
  
  
  # ----Perform Gauss-Seidel Iterations for NB > 2----
  Mlist <- NA 
  cat("\n\n")
  
  if( NB > 2 ){
        
     # ----____Create Initial Mlist, Cij, Bij, ---- 
     # needed for Chi square when NB = 2 and 
     # Gauss Seidel when NB > 2
     Mlist <- fncMlist(NB, NVB, mGamma, NFac) 

     # ----____Create BlockIJ----
     BlockIJ <- fncBlockIJ(NB, NVB)  
          
     #Estimate LambdaHat 
     LambdaHat <- fncGaussSeidel(Tmat, NVB, NB, V, Mlist, 
                                 NFac, BlockIJ, epsilon, PrintLevel)
    
     # Update Mlist a final time for computing X^2 values
     # mGamma is global from fncUpdateGamma function
     Mlist <- fncMlist(NB, NVB, mGamma, NFac) 
   } ##END if NB > 2
  


  #----Rotate Local Solutions----         
   ListRotatedSolutions <- fncRotateLocalSolutions(LambdaHat, cnRotate) 
  #----Order Local Solutions----
   OutfncOrderUniqueSolutions <-fncOrderUniqueSolutions(ListRotatedSolutions,cnRotate) 
      uniqueSolutions <-OutfncOrderUniqueSolutions$uniqueSolutions
      nGrp <-OutfncOrderUniqueSolutions$nGrp
      localGrps <- OutfncOrderUniqueSolutions$localGrps
      UnSpunPosition <- OutfncOrderUniqueSolutions$UnSpunPosition
  
     # The global solution: (smallest rotation complexity value 
     # is always in position [[1]])
     fout <- uniqueSolutions[[1]]
  
  
  # ----Compute Model Fit ----
    FitOut <- fncModelFit(R, 
                          NB,
                          NVB,
                          Mlist,
                          V,
                          numFactors, 
                          Lambda = fout$loadings, 
                          Phi = fout$Phi,
                          SampleSize)
    
    fitVec <- list(SampleSize = SampleSize,
                   ChiSq = FitOut$ChiSq,
                   DF = FitOut$DF,
                   pvalue = FitOut$pvalue,
                   AIC = FitOut$AIC,
                   BIC = FitOut$BIC,
                   RMSEA = FitOut$RMSEA,
                   rhoSq = FitOut$rhoSq)
    
    Rhat <-FitOut$Rhat
    Resid <-FitOut$Resid
    
  
 #----OUTPUT ----  
  

  Output <- list(loadings              = fout$loadings, # this is min complexity solutution
                 Phi                   = fout$Phi,
                 fit                   = fitVec,
                 R                     = R,
                 Rhat                  = Rhat,
                 Resid                 = Resid,
                 facIndeterminacy      = fout$facIndeter,
               #  loadingsSE            = loadSE,
               #  CILevel               = CILevel,
               #  loadingsCIupper       = loadCI.upper,
               #  loadingsCIlower       = loadCI.lower,
               #  PhiSE                 = phiSE,
               #  PhiCIupper            = phiCI.upper,
               #  PhiCIlower            = phiCI.lower,
               #  facIndeterminacySE    = FISE,
                  localSolutions        = uniqueSolutions,
                  numLocalSets          = nGrp,
                  localSolutionSets     = localGrps,
                  rotate                = rotate,
               #  loadingsArray         = loadArray,
               #  PhiArray              = phiArray,
               #  facIndeterminacyArray = FIArray,
                  rotateControl         = cnRotate,
                  unSpunSolution        = uniqueSolutions[[UnSpunPosition]],
                  Call = Call) #END list( . . .
  
    ## ----____Add Names and Class to the Output----
    ## Add labels to the resulting output  
    fout <- fncAddDimNames(Output)
  
  ## Function return:
  class(fout) <- 'faMB'
  fout
} #END faMB
