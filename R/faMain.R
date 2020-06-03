#' Automatic Factor Rotation from Random Configurations with Bootstrap Standard Errors
#' 
#' This function conducts factor rotations (using the \pkg{GPArotation} package) 
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
#' @param numFactors (Numeric) The number of factors to extract for subsequent rotation.
#' @param facMethod (Character) The method used for factor extraction 
#' (\code{\link{faX}}). The supported options are "fals" for unweighted least 
#' squares, "faml" for maximum likelihood, "fapa" for iterated principal axis 
#' factoring, "faregLS" for regularized least squares,
#' "faregML" for regularized maximum likelihood, and "pca" for principal components 
#'  analysis. The default method  is "fals". 
#' \itemize{
#'   \item \strong{"fals"}: Factors are extracted using the unweighted least 
#'   squares estimation procedure using the \code{\link{fals}} function.
#'   \item \strong{"faml"}: Factors are extracted using the maximum likelihood 
#'   estimation procedure using the \code{\link[stats]{factanal}} function.
#'   \item \strong{"fapa"}: Factors are extracted using the iterated principal 
#'   axis factoring estimation procedure using the \code{\link{fapa}} function.
#'   \item \strong{"faregLS"}: Factors are extracted using regularized 
#'   least squares factor analysis using the \code{\link{fareg}} function. 
#'   \item \strong{"faregML"}: Factors are extracted using regularized 
#'   maximum likelihood factor using the \code{\link{fareg}} function. 
#'   \item \strong{"pca"}: Principal components are extracted. 
#' }
#' @param urLoadings (Matrix) An unrotated factor-structure matrix to be rotated.
#' @param rotate (Character) Designate which rotation algorithm to apply. The 
#' following are available rotation options: "oblimin", "quartimin", "targetT", 
#' "targetQ", "oblimax", "entropy", "quartimax", "varimax", "simplimax", 
#' "bentlerT", "bentlerQ", "tandemI", "tandemII", "geominT", "geominQ", "cfT", 
#' "cfQ", "infomaxT", "infomaxQ", "mccammon", "bifactorT", "bifactorQ", and 
#' "none". Defaults to rotate = "oblimin". See \pkg{GPArotation} package for more 
#' details. Note that rotations ending in "T" and "Q" represent orthogonal and 
#' oblique rotations, respectively.
#' @param targetMatrix (Matrix) This argument serves two functions. First, if a 
#' user has requested either a "targetT" or "targetQ' rotation, then 
#'  the target matrix is used to conduct a fully or partially
#' specified target rotation. In the latter case,  freely estimated factor 
#' loadings are designated by "NA" values and rotation will be conducted using  
#' Browne's (1972a, 1972b, 2001) method for a partially-specified 
#' target rotation. Second, if any other rotation option is chosen then all 
#' rotated loadings matrices (and assorted output) will be aligned 
#' (but not rotated) with the target solution. 
#' @param bootstrapSE (Logical) Computes bootstrap standard errors. All bootstrap 
#' samples are aligned to the global minimum solution. Defaults to 
#' bootstrapSE = FALSE (no standard errors). 
#' @param numBoot (Numeric) The number bootstraps. Defaults to numBoot = 1000.
#' @param CILevel (Numeric) The confidence level (between 0 and 1) of the bootstrap 
#' confidence interval. Defaults to CILevel = .95.
#' @param Seed (Numeric) Starting seed for reproducible bootstrap results and factor rotations. 
#' Defaults to Seed = 1.
#' @param digits (Numeric) Rounds the values to the specified number of decimal 
#' places. Defaults to digits = NULL (no rounding).
#' @param faControl (List) A list of optional parameters passed to the factor 
#' extraction (\code{\link{faX}}) function.
#' \itemize{
#'   \item \strong{treatHeywood}: (Logical) In \code{fals}, if treatHeywood is 
#'   true, a penalized least squares function is used to bound the communality 
#'   estimates below 1.0. Defaults to treatHeywood = TRUE.
#'   \item \strong{nStart}: (Numeric) The number of starting values to be tried 
#'   in \code{faml}. Defaults to nStart = 10.
#'   \item \strong{start}: (Matrix) NULL or a matrix of starting values, each column 
#'   giving an initial set of uniquenesses. Defaults to start = NULL. 
#'   \item \strong{maxCommunality}: (Numeric) In \code{faml}, set the maximum 
#'   communality value for the estimated solution. Defaults to maxCommunality = .995.
#'   \item \strong{epsilon}: (Numeric) In \code{fapa}, the numeric threshold 
#'   designating when the algorithm has converged. Defaults to epsilon = 1e-4.
#'   \item \strong{communality}: (Character) The method used to estimate the 
#'   initial communality values in \code{fapa}. Defaults to communality = 'SMC'.
#'   \itemize{
#'     \item \strong{"SMC"}: Initial communalities are estimated by taking the 
#'     squared multiple correlations of each indicator after regressing the 
#'     indicator on the remaining variables.
#'     \item \strong{"maxr"}: Initial communalities equal the largest 
#'     (absolute value) correlation in each column of the correlation matrix.
#'     \item \strong{"unity"}: Initial communalities equal 1.0 for all variables.
#'   }
#'   \item \strong{maxItr}: (Numeric) In \code{fapa}, the maximum number of 
#'   iterations to reach convergence. Defaults to maxItr = 15,000.
#' }
#' 
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
#' 
#' @param ... Values to be passed to the \code{\link[stats]{cor}} function.
#' \itemize{
#'   \item \strong{use}: (Character) A character string giving a method for 
#'   computing correlations in the presence of missing values: "everything" 
#'   (the default), "all.obs", "complete.obs", "na.or.complete", or 
#'   "pairwise.complete.obs".
#'   \item \strong{method}: (Character) A character string indicating which 
#'   correlation coefficient is to be computed: "pearson" (the default), 
#'   "kendall", or "spearman". 
#'   \item \strong{na.rm}: (Logical) Should missing values be removed (TRUE) 
#'   or not (FALSE)?
#' }
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
#' @return The \code{faMain} function will produce a lot of output in addition 
#' to the rotated factor pattern matrix and the factor correlations.
#' \itemize{
#'   \item \strong{R}: (Matrix) Returns the correlation matrix, useful when raw data are supplied.
#'   \item \strong{loadings}: (Matrix) The rotated factor solution with the 
#'   lowest evaluated discrepancy function. This solution has the lowest 
#'   discrepancy function \emph{of the examined random starting configurations}. 
#'   It is not guaranteed to find the "true" global minimum. Note that multiple
#'    (or even all) local solutions can have the same discrepancy functions.
#'   \item \strong{Phi}: (Matrix) The factor correlations of the rotated factor 
#'   solution with the lowest evaluated discrepancy function (see Details).
#'   \item \strong{facIndeterminacy}: (Vector) A vector (with length equal to the number of factors)
#'   containing Guttman's (1955) index of factor indeterminacy for each factor. 
#'   \item \strong{h2}: (Vector) The vector of final communality estimates. 
#'   \item \strong{loadingsSE}: (Matrix) The matrix of factor-loading standard 
#'   errors across the bootstrapped factor solutions. Each matrix element is 
#'   the standard deviation of all bootstrapped factor loadings for that element position.
#'   \item \strong{CILevel} (Numeric) The user-defined confidence level (between 0 and 1) of the bootstrap 
#'    confidence interval. Defaults to CILevel = .95.
#'   \item \strong{loadingsCIupper}: (Matrix) Contains the upper confidence 
#'   interval of the bootstrapped factor loadings matrix. The confidence 
#'   interval width is specified by the user.
#'   \item \strong{loadingsCIlower}: (Matrix) Contains the lower confidence 
#'   interval of the bootstrapped factor loadings matrix. The confidence 
#'   interval width is specified by the user.
#'   \item \strong{PhiSE}: (Matrix) The matrix of factor correlation standard 
#'   errors across the bootstrapped factor solutions. Each matrix element is 
#'   the standard deviation of all bootstrapped factor correlations for that element position.
#'   \item \strong{PhiCIupper}: (Matrix) Contains the upper confidence interval 
#'   of the bootstrapped factor correlation matrix. The confidence interval 
#'   width is specified by the user.
#'   \item \strong{PhiCIlower}: (Matrix) Contains the lower confidence interval 
#'   of the bootstrapped factor correlation matrix. The confidence interval 
#'   width is specified by the user.
#'   \item \strong{facIndeterminacySE}: (Matrix) A row vector containing the 
#'   standard errors of Guttman's (1955) factor indeterminacy indices across the 
#'   bootstrap factor solutions. 
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
#'   \item \strong{loadingsArray}: (Array) Contains an array of all bootstrapped 
#'   factor loadings. The dimensions are factor indicators, factors, and the 
#'   number of bootstrapped samples (representing the row, column, and depth, respectively).
#'   \item \strong{PhiArray}: (Array) Contains an array of all bootstrapped 
#'   factor correlations. The dimension are the number of factors, the number 
#'   of factors, and the number of bootstrapped samples (representing the row,
#'    column, and depth, respectively).
#'   \item \strong{facIndeterminacyArray}: (Array) Contains an array of all 
#'   bootstrap factor indeterminacy indices. The dimensions are 1, the number 
#'   of factors, and the number of bootstrap samples (representing the row, 
#'   column, and depth order, respectively).
#'   \item \strong{faControl}: (List) A list of the control parameters passed 
#'   to the factor extraction (\code{\link{faX}}) function.
#'   \item \strong{faFit}: (List) A list of additional output from the factor
#'   extraction routines. 
#'  \itemize{
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
#'     \item \strong{convergedX}: (Logical) TRUE if the factor 
#'     \strong{extraction} routine converged. 
#'     \item \strong{convergedR}: (Logical) TRUE if the factor \strong{rotation} 
#'     routine converged (for the local solution with the minimum discrepancy 
#'     value).
#'   }
#'   \item \strong{rotateControl}: (List) A list of the control parameters 
#'   passed to the rotation algorithm.
#'   \item \strong{unSpunSolution}: (List) A list of output parameters (e.g., loadings, Phi, etc) from 
#'   the rotated solution that was obtained by rotating directly from the unrotated (i.e., unspun) common factor orientation. 
#'   \item \strong{targetMatrix} (Matrix) The input target matrix if supplied by the user.
#'   \item \strong{Call}: (call) A copy of the function call.
#' }
#' 
#' @references Browne, M. W.  (1972).  Oblique rotation to a partially specified 
#' target.  \emph{British Journal of Mathematical and Statistical Psychology, 25},(1), 
#' 207-212.  
#' @references Browne, M. W. (1972b). Orthogonal rotation to a partially specifed
#' target. \emph{British Journal of Statistical Psychology, 25},(1), 115-120.
#' @references Browne, M. W. (2001). An overview of analytic rotation in 
#' exploratory factor analysis. \emph{Multivariate Behavioral Research, 36}(1), 111-150.
#' @references Cureton, E. E., & Mulaik, S. A. (1975). The weighted varimax 
#' rotation and the promax rotation. \emph{Psychometrika, 40}(2), 183-195.
#' @references Guttman, L. (1955). The determinacy of factor score matrices with 
#' implications for five other basic problems of common factor theory. 
#' \emph{British Journal of Statistical Psychology, 8}(2), 65-81.
#' @references Jung, S. & Takane, Y.  (2008).  Regularized common factor analysis.  
#' \emph{New Trends in Psychometrics}, 141-149. 
#' @references Mansolf, M., & Reise, S. P. (2016). Exploratory bifactor 
#' analysis: The Schmid-Leiman orthogonalization and Jennrich-Bentler 
#' analytic rotations. \emph{Multivariate Behavioral Research, 51}(5), 698-717.
#' @references Rozeboom, W. W. (1992). The glory of suboptimal factor rotation: 
#' Why local minima in analytic optimization of simple structure are more 
#' blessing than curse. \emph{Multivariate Behavioral Research, 27}(4), 585-599.
#' @references Zhang, G. (2014). Estimating standard errors in exploratory factor 
#' analysis. \emph{Multivariate Behavioral Research, 49}(4), 339-353.
#'
#' @family Factor Analysis Routines
#'
#' @author
#' \itemize{
#'   \item Niels G. Waller (nwaller@umn.edu)
#'   \item Casey Giordano (Giord023@umn.edu)
#'   \item The authors thank Allie Cooperman and Hoang
#'    Nguyen for their help implementing the standard error estimation and the 
#'    Cureton-Mulaik standardization procedure.
#'}
#'
#' @examples
#' ## Example 1
#' 
#' ## Generate an orthgonal factor model
#' lambda <- matrix(c(.41, .00, .00,
#'                    .45, .00, .00,
#'                    .53, .00, .00,
#'                    .00, .66, .00,
#'                    .00, .38, .00,
#'                    .00, .66, .00,
#'                    .00, .00, .68,
#'                    .00, .00, .56,
#'                    .00, .00, .55),
#'                  nrow = 9, ncol = 3, byrow = TRUE)
#'
#' ## Generate factor correlation matrix
#' Phi <- matrix(.50, nrow = 3, ncol = 3)
#' diag(Phi) <- 1
#'
#' ## Model-implied correlation matrix
#' R <- lambda %*% Phi %*% t(lambda)
#' diag(R) <- 1
#'
#' ## Load the MASS package to create multivariate normal data
#' library(MASS)
#' 
#' ## Generate raw data to perfectly reproduce R
#' X <- mvrnorm(Sigma = R, mu = rep(0, nrow(R)), empirical = TRUE, n = 300)
#' 
#'\dontrun{
#' ## Execute 50 promax rotations from a least squares factor extraction
#' ## Compute 100 bootstrap samples to compute standard errors and 
#' ## 80 percent confidence intervals
#' Out1 <- faMain(X             = X,
#'                numFactors    = 3,
#'                facMethod     = "fals",
#'                rotate        = "promaxQ",
#'                bootstrapSE   = TRUE,
#'                numBoot       = 100,
#'                CILevel       = .80,
#'                faControl     = list(treatHeywood = TRUE),
#'                rotateControl = list(numberStarts = 2,  
#'                                     power        = 4,
#'                                     standardize  = "Kaiser"),
#'                digits        = 2)
#' Out1[c("loadings", "Phi")] 
#'}
#'
#' ## Example 2
#'
#' ## Load Thurstone's (in)famous box data
#' data(Thurstone, package = "GPArotation")
#' 
#' ## Execute 5 oblimin rotations with Cureton-Mulaik standardization 
#' Out2 <- faMain(urLoadings    = box26,
#'                rotate        = "oblimin",
#'                bootstrapSE   = FALSE,
#'                rotateControl = list(numberStarts = 5,
#'                                     standardize  = "CM",
#'                                     gamma        = 0,
#'                                     epsilon      = 1e-6),
#'                digits        = 2)
#'                
#' Out2[c("loadings", "Phi")]     
#' 
#' ## Example 3
#' 
#' ## Factor matrix from Browne 1972
#' lambda <- matrix(c(.664,  .322, -.075,
#'                    .688,  .248,  .192,
#'                    .492,  .304,  .224,
#'                    .837, -.291,  .037,
#'                    .705, -.314,  .155,
#'                    .820, -.377, -.104,
#'                    .661,  .397,  .077,
#'                    .457,  .294, -.488,
#'                    .765,  .428,  .009), 
#'                  nrow = 9, ncol = 3, byrow = TRUE)   
#'                  
#' ## Create partially-specified target matrix
#' Targ <- matrix(c(NA, 0,  NA,
#'                  NA, 0,  0,
#'                  NA, 0,  0,
#'                  NA, NA, NA,
#'                  NA, NA, 0,
#'                  NA, NA, NA,
#'                  .7, NA, NA,
#'                  0,  NA, NA,
#'                  .7, NA, NA), 
#'                nrow = 9, ncol = 3, byrow = TRUE)  
#'                
#' ## Perform target rotation              
#' Out3 <- faMain(urLoadings   = lambda,
#'                rotate       = "targetT",
#'                targetMatrix = Targ,
#'                digits       = 3)$loadings
#' Out3
#' @import GPArotation
#' @import stats
#' @export

faMain <-
  function(X             = NULL,      ## Raw data matrix or data.frame
           R             = NULL,      ## Correlation matrix for extraction
           n             = NULL,      ## Sample Size for faml
           numFactors    = NULL,      ## Number of factors to extract
           facMethod     = "fals",    ## Factor extraction method
           urLoadings    = NULL,      ## Unrotated loading matrix
           rotate        = "oblimin", ## Which rotation to use
           targetMatrix  = NULL,      ## If targetT or targetQ, matrix of specified loadings
           bootstrapSE   = FALSE,     ## Whether to bootstrap stand. errors
           numBoot       = 1000,      ## Number of bootstrapped samples drawn
           CILevel       = .95,       ## Bootstrap SE confidence level
           Seed          = 1,         ## Seed for bootstrap reproducibility
           digits        = NULL,      ## Round all output
           faControl     = NULL,      ## Control factor extraction parameters
           rotateControl = NULL,      ## Control rotation tuning parameters
           ...) {
    
    ## ~~~~~~~~~~~~~~~~~~ ##
    #### Error checking ####
    ## ~~~~~~~~~~~~~~~~~~ ##
    
    if( !is.null(R) ) {
       # Check if only row or col names given
       # then remove names
       rowOrCol <- which.max(c(length(rownames(R)),
                            length(colnames(R))))
       Vnames <- dimnames(R)[[rowOrCol]]
       #Put row and col names back on matrix
       dimnames(R) <- list(Vnames,Vnames)
    }# END  if(!is.null(R))   
    
    
    
    FLAGurLoadings <- FALSE
    if(!is.null(urLoadings))  FLAGurLoadings <- TRUE
    
    ## Check rotate
    
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
    
    ## Create nVar object, used in defining 'k' param in cnRotate list
    if (!is.null(X))          nVar <- ncol(X)
    if (!is.null(R))          nVar <- ncol(R)
    if (!is.null(urLoadings)) nVar <- nrow(urLoadings)
    
    ## Check R
    
    ## Make sure either corr mat or factor structure is supplied
    if ( is.null(R) & is.null(urLoadings) & is.null(X) ) {
      stop("Must specify 'X', 'R', or 'urLoadings' arguments.")
    } # END if ( is.null(R) & is.null(urLoadings) & is.null(X) )
    
    ## If factor analysis is run, must specify numFactors
    if ( is.null(urLoadings) & is.null(numFactors) ) {
      stop("The user needs to specify the number of factors to extract.")
    } # END if ( is.null(urLoadings) & is.null(numFactors) )
    
    ## Check X (if specified)
    
    ## Only check the following if X is specified
    if ( !is.null(X) ) {
      
      
      ## Missingness in data?
      if (nrow(X) != nrow(X[complete.cases(X),])) {
        warning("There are missing values in the data frame. See details for how missing data is handled in the cor function.")
      } # END if (nrow(X) != nrow(X[complete.cases(X),]))
      
      ## Class must be matrix or DF to analyze via cor() function
      if ( !any( class(X) %in% c("matrix", "data.frame", "loadings") ) ){
        stop("'X' must be of class matrix, data.frame, or loadings.")
      } # END if ( class(X) %in% c("matrix", "data.frame", "loadings"))
      
    } # END if ( !is.null(X) )
    
    ## Digits
    
    ## If a digits argument is not supplied, use R's default option
    if ( is.null(digits) ) {
      digits <- options()$digits
    } # END if ( is.null(digits) )
    
    ## BootstrapSE
    
    ## Must give the raw data to perform bootstraps
    if ( bootstrapSE == TRUE && is.null(X)) {
      stop("The raw data are required to compute the bootstrapped standard errors.")
    } # END if ( bootstrapSE == TRUE && is.null(X))
    
    ## rotateControl
    
    #### ------- DEFINE cnRotate -------- ####
    
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
      
    } # END if ( is.null(rotateControl) )
    
    ## Change the default values based on user-specified rotateControl arguments
    cnRotate[names(rotateControl)] <- rotateControl
    
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    #### Define Utility Function(s) ####
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    
    ## Generate a random orthonormal starting matrix
    randStart <- function(dimension) {
      qr.Q(qr(matrix(rnorm(dimension^2), dimension, dimension )))
    } # END randStart <- function(dimension) {
    
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
             "targetT" = {
               GPArotation::targetT(lambda,
                                    Tmat      = spinMatrix,
                                    Target    = targetMatrix,
                                    normalize = FALSE,
                                    eps       = cnRotate$epsilon,
                                    maxit     = cnRotate$maxItr)
             },
             "targetQ" = {
               GPArotation::targetQ(lambda,
                                    Tmat      = spinMatrix,
                                    Target    = targetMatrix,
                                    normalize = FALSE,
                                    eps       = cnRotate$epsilon,
                                    maxit     = cnRotate$maxItr)
             },
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
      
      ## Try to invert, if not invertible, gives class 'try-error
      Rinv  <-try(solve(SampCorr), silent = TRUE)
      ## If non-invertible, return NAs instead of returning an error
      if ( any( class(Rinv) %in% "try-error") ) {
        # warning("\n\nEncountered a singular R matrix when computing factor indeterminancy values\n")
        return( rep(NA, ncol(Lambda) ) )
      } # END if ( any( class(Rinv) %in% "try-error") ) 
      ## Factor indeterminacy solution
      FI <- sqrt( diag( t(facStruct) %*% Rinv %*% facStruct))
      if(max(FI)>1) FI <- rep(NA, length(FI))
      FI
    } # END GuttmanIndices
    
    #### ------- VARIABLE NAME RETENTION -------- ####
    
    ## Extract names of the indicators to rename final output
    varNames <- NULL 
    
    ## If variable names are provided, retain them
    if ( !is.null(X) )          varNames <- colnames(X)
    if ( !is.null(R) )          varNames <- colnames(R)
    if ( !is.null(urLoadings) ) varNames <- rownames(urLoadings)
    
    ## If variable names is NULL, assign them as var 1 to var nItems
    if ( is.null(varNames) ) {
      
      ## Specify the number of indicators based on the function inputs
      if ( !is.null(X) )         nItem <- ncol(X)
      if ( !is.null(R) )         nItem <- ncol(R)
      if ( !is.null(urLoadings)) nItem <- nrow(urLoadings)
      
      varNames <- paste0("var", seq_len(nItem))
    } # END if ( is.null(varNames) ) 
    
    #### ------- BEGIN FUNCTION -------- ####
    
    ## Return a copy of the call function. 
    CALL <- match.call()
    
    ## Compute the correlation matrix
    if ( !is.null(X) ) {
      
      ## Create correlation matrix (default)
      R <- cor(X, ...)
      
    } # END if ( is.null(R) )
    
    #### ------- FACTOR EXTRACTION -------- ####
    
    ## If unrotated loadings specified, 
    if ( !is.null(urLoadings) ) {
      
      ## September 3, 2019  Check Changes
      # needed if no factor extraction occurs
      faOut <- vector(mode="list")
      faOut$faControl <- NULL
      
      ## No model fit stuff
      faModelFit <- NULL
      
      ## Compute communalities before rotation
      faXh2 <- apply(urLoadings^2, 1, sum)
      
    } # END if ( !is.null(urLoadings) ) 
    
    ## If urLoadings is not supplied, extract it
    if ( is.null(urLoadings) ) {
      
      ## Call faX function to extract factor structure matrix
      ## September 3, 2019  Check Changes
      faOut <- vector(mode="list")
      faOut$faControl <- NULL
      
      faOut <- faX(R          = R,
                   n          = n,
                   numFactors = numFactors,
                   facMethod  = facMethod,
                   faControl  = faControl)
      
      
      
      ## Extract factor loadings
      urLoadings <- faOut$loadings[]
      
      ## Communalities from factor extraction
      faXh2 <- faOut$h2
      
      ## Extract all other output
      
      convergedPos <- which(names(faOut$faFit) %in% "converged")
      names(faOut$faFit)[convergedPos] <- "convergedX"
      faModelFit <- faOut$faFit
      
    } # END if ( is.null(urLoadings) )
    
    ## If urLoadings is provided and numFactors isn't, specify numFactors
    if ( is.null(numFactors) && !is.null(urLoadings) ) {
      numFactors <- ncol(urLoadings)
    } # END if ( is.null(numFactors) && !is.null(urLoadings) )
    
    #### ------- STANDARDIZATION -------- ####
    
    ## Call standardize function
    Stnd <- faStandardize(method = cnRotate$standardize,
                          lambda = urLoadings)
    
    ## Extract DvInv for later unstandardization
    lambda <- Stnd$lambda
    
    #### ------- ROTATION -------- ####
    
    
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
    starts <- lapply(starts, function(x) {
      
      ## Rotate the standardized factor structure matrix
      rotatedLambda <- internalRotate(lambda        = lambda,
                                      rotation      = rotate,
                                      spinMatrix    = x,
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
    
    #### ------- ORDER UNIQUE SOLUTIONS -------- ####
    
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
      
      #### ------- FACTOR SORT -------- ####
      
      
      ## June 27, 2019: Do not  order factors when Procrustes rotation chosen
      ##
      ## sort loadings and Phi matrix for non Procrstes
      if( rotate != "targetT" & rotate != "targetQ"){
        sortedSols <- 
          orderFactors(Lambda  = selected$loadings,
                       PhiMat  = selected$Phi,
                       salient = .01,
                       reflect = TRUE)
        
        ## Overwrite "selected" list
        selected$loadings <- sortedSols$Lambda
        selected$Phi      <- sortedSols$PhiMat
      }
      
      
      
      
      ## Select the factor loadings
      uniqueSolutions[[iternumStart]]$loadings <- selected$loadings
      
      ## Select factor correlations
      uniqueSolutions[[iternumStart]]$Phi <- selected$Phi
      
      ## Select complexity values
      uniqueSolutions[[iternumStart]]$RotationComplexityValue <- ComplexityFunc[num]
      
      #### ----- FACTOR INDETERMINACY ----- ####
      ## if the user provides urloadings then do not 
      ## compute factor indeterminancy values
    
      if( !isTRUE(FLAGurLoadings) ){
        ## Guttman's factor indeterminacy indices
        uniqueSolutions[[iternumStart]]$facIndeterminacy <- 
          GuttmanIndices(Lambda = selected$loadings, 
                         PhiMat = selected$Phi,
                         SampCorr = R ) 
      }else{
        uniqueSolutions[[iternumStart]]$facIndeterminacy <- rep(NA, numFactors)
      }
      
      # END if(!is.null((R)))
      
      ## Did the local optima solution converge
      uniqueSolutions[[iternumStart]]$RotationConverged <- selected$convergence
      
    } # END for (iternumStart in 1:length(sortedDiscOrder))
    
    ## Extracted rotation complexity values from uniqueSolutions to find solution sets
    DisVal <- unlist(lapply(uniqueSolutions, function(x) x$RotationComplexityValue))
    
    ## Round the discrepancy values to specified number of digits
    DisVal <- round(x      = DisVal,
                    digits = nchar(cnRotate$epsilon))
    
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
    facIndeter <- data.frame("FI" = uniqueSolutions[[1]]$facIndeterminacy)

    ## If factor indeter. not computed, give NA values
    if ( is.null(facIndeter) ) facIndeter <- rep(NA, ncol(minLambda))
    
    ## If target matrix given then align final factor solution 
    ## to the target
    
    if ( !is.null(targetMatrix) && 
         rotate %in% c("targetT", "targetQ", "none") == FALSE ) {
      minAlign   <- faAlign(F1          = targetMatrix, 
                            F2          = minLambda, 
                            Phi2        = minPhi,
                            MatchMethod = "LS")
      minLambda  <- minAlign$F2
      minPhi     <- minAlign$Phi2
      facIndeter <- data.frame("FI" = facIndeter[ minAlign$FactorMap[2, ], ])
    }
    
    
    
    ## CG EDITS (30 sept 19): Moved communality computation to before bootstraps
    ## Compute communalities. If orthogonal model, Phi is identity matrix
    communalityDF <- 
      data.frame("h2" = diag(minLambda %*% minPhi %*% t(minLambda)))
    
    #### -------- BOOTSTRAP SETUP -------- ####
    
    ## If true, compute bootstrap standard errors
    if (bootstrapSE == TRUE) {
      
      # Number of subjects
      nSubj <- nrow(X)
      
      # Lists/variables to hold output
      ## Hold factor loadings from bootstraps
      fList <-
        ## Hold factor correlations from bootstraps
        phiList <- 
        ## Hold factor indeterminacy indices from bootstraps
        FIList <- 
        ## Hold communality estimates from bootstraps
        h2List <- vector("list", numBoot)
      
      ## Number of rows, used below for bootstrap sampling
      rows <- 1:nSubj
      
      #### ----____BOOTSTRAP FOR LOOP -------- ####
      
      ## Analyses on 'numBoot' number of random resamples of X matrix
      for (iSample in seq_len(numBoot)) {
        
        ## Set the seed for reproducibility
        set.seed(iSample + Seed)
        
        ## Resample (with replacement) from X raw data matrix
        bsSample <- sample(x       = rows,
                           size    = nSubj,
                           replace = TRUE)
        
        ## Create correlation matrix from the bootstrap sample
        # Rsamp <- cor(X[bsSample, ], ...)
        # if missing data present in the original sample, use
        # missing data method passed to faMain fnc 
        
        Rsamp <- cor(X[bsSample, ], ...)

        
        ## Extract unrotated factors using resampled data matrix
        bsLambda <- faX(R          = Rsamp,
                        numFactors = numFactors,
                        facMethod  = facMethod,
                        faControl  = faControl)$loadings
        
        ## Conduct standardization
        bsStnd <- faStandardize(method = cnRotate$standardize,
                                lambda = bsLambda)
        
        ## Extract the standardized bootstrapped (unrotated) factor loadings
        bsLambda <- bsStnd$lambda
        
        ## Find the "global" min out of all the random starting configurations
        bsStarts <- vector("list", cnRotate$numberStarts)
        
        bsStarts <- lapply(bsStarts, function(x) randStart(dimension = numFactors))
        
        ## Conduct rotations from random start values
        bsStarts <- lapply(bsStarts, function(x) {
          
          ## Rotate the bootstrapped samples
          bsRotated <- internalRotate(lambda        = bsLambda,
                                      rotation      = rotate,
                                      spinMatrix    = x,
                                      rotateControl = cnRotate)
          
          ## If an orthogonal model, turn Phi into an identity matrix
          if ( is.null(bsRotated$Phi) ) {
            
            ## Identity matrix
            bsRotated$Phi <- diag(numFactors)
            
          } # END if ( is.null(bsRotated$Phi) ) 
          
          ## Return all output, not just Phi matrix
          bsRotated
          
        }) # END bsStarts <- lapply(bsStarts, function(x)
        
        ## Find minimum of bsSolutions
        ## Evaluate the minimum disc functions to find smallest value
        if ( rotate != "none" ) {
          
          ## For each random start, find the evaluated discrepancy function
          bootstrapComplexityFunc <- sapply(bsStarts, function(attempt) min(attempt$Table[, 2]))
          
        }# END if ( rotate != "none" ) 
        
        ## If no rotation, no discrepancy value to evaluate
        if (rotate == "none") bootstrapComplexityFunc <- rep(1, cnRotate$numberStarts)
        
        ## Of all random configs, determinine which has the lowest criterion value
        bsLambda <- bsStarts[[which.min(bootstrapComplexityFunc)]]$loadings
        bsPhi    <- bsStarts[[which.min(bootstrapComplexityFunc)]]$Phi
        
        ## Unstandardize the minimum rotated solution for this bootstrap sample
        bsLambda <- bsStnd$DvInv %*% bsLambda
        
        ## Pre-allocate list
        Aligned <- list()
        
        ## Factor alignment does not work on 1 factor (and doesn't make sense)
        if (numFactors > 1) {
          ## Align bootstrap sample with global minimum found earlier
          Aligned <- faAlign(F1   = minLambda,
                             F2   = bsLambda,
                             Phi2 = bsPhi)
          
          ## Extract aligned elements
          AlignedLambda <- Aligned$F2
          AlignedPhi    <- Aligned$Phi2
        } # END if (numFactors > 1) 
        
        ## Cannot align 1-factor model, but can reflect the factor
        if (numFactors == 1) {
          if ( sum(bsLambda) < 0 ) bsLambda <- bsLambda * -1
          
          ## Define newly aligned elements
          AlignedLambda <- bsLambda
          AlignedPhi    <- bsPhi
        } # END if (numFactors == 1)
        
        ## Save loadings as the bootstrapped factor loadings
        fList[[iSample]]   <- AlignedLambda
        
        ## Save correlations as the bootstrapped factor correlations
        phiList[[iSample]] <- AlignedPhi
        
        ## Save factor indeterminacy indices
        FIList[[iSample]] <- GuttmanIndices(Lambda = AlignedLambda,
                                            PhiMat = AlignedPhi,
                                            SampCorr = Rsamp)
        
        ## CG EDITS (30 sept 19): Added list of computed h^2 values
        ## Save communality estimates
        h2List[[iSample]] <- diag(AlignedLambda %*% AlignedPhi %*% t(AlignedLambda))
        
      } # END for (iSample in seq_len(numBoot))
      
      #### -----____ BootStrap STANDARD ERRORS ----- ####
      
      # Convert list of matrices into array
      loadArray <- array(unlist(fList), c(nVar, numFactors, numBoot))
      phiArray  <- array(unlist(phiList), c(numFactors, numFactors, numBoot))
      FIArray   <- array(unlist(FIList), c(1, numFactors, numBoot))
      ## CG EDITS (30 sept 19): Added array for h2 values
      h2Array   <- array(unlist(h2List), c(nVar, 1, numBoot))
      
      
      # Bootstrap standard errors for factor loadings
      loadSE <- apply(loadArray, 1:2, sd)
      loadSE <- round(loadSE, digits)
      
      # Bootstrap standard errors for factor correlations
      phiSE <- apply(phiArray, 1:2, sd)
      phiSE <- round(phiSE, digits)
      
      
      ## CG EDITS (30 sept 19): Added bootstrap SEs for communalities
      h2SE <- apply(h2Array, 1, sd)
      
      
      ## ----- CONFIDENCE INTERVALS ----- ##
      
      # Set alpha level
      alpha <- 1 - CILevel
      
      # Confidence interval matrices for factor loadings
      loadCI.upper <- apply(loadArray, 1:2, quantile, 1 - (alpha / 2))
      loadCI.lower <- apply(loadArray, 1:2, quantile, (alpha / 2))
      
      ## Round loadings
      loadCI.upper <- round(loadCI.upper, digits)
      loadCI.lower <- round(loadCI.lower, digits)
      
      # Confidence interval matrices for factor correlations
      phiCI.upper <- apply(phiArray, 1:2, quantile, 1 - (alpha / 2))
      phiCI.lower <- apply(phiArray, 1:2, quantile, (alpha / 2))
      
      ## Round Phi
      phiCI.upper <- round(phiCI.upper, digits)
      phiCI.lower <- round(phiCI.lower, digits)
      
      ## CG EDITS (30 sept 19): Added CIs for fact indeterminacy estimates
      ## ----____Bootstrap SE FI----
      if(!any(is.na(FIArray))){
          FISE <- apply(FIArray, 1:2, sd)
         # FISE <- round(FISE, digits)
          FICI.upper <- apply(FIArray, 2, quantile, 1 - (alpha / 2))
          FICI.lower <- apply(FIArray, 2, quantile, (alpha / 2))
      
          ## CG EDITS (30 sept 19): Create data frame of fac indeterminacy estimates
          facIndeter <- data.frame(facIndeter,
                               t(FISE),
                               FICI.lower,
                               FICI.upper)
          ## CG EDITS (30 sept 19): Add column names to data frame of FI estimates
          colnames(facIndeter) <- c("FI",
                                    "SE",
                                    paste0((alpha / 2) * 100, "th percentile"),
                                    paste0((1 - (alpha / 2)) * 100, "th percentile")) 
          rownames(facIndeter) <-  paste0("f", 1:numFactors)
      }else{
        facIndeter <- NA
        FISE <- NA
      }    
      
    
      
      ## CG EDITS (30 sept 19): Added CIs for communality estimates
      h2CI.upper <- apply(h2Array, 1, quantile, 1 - (alpha / 2))
      h2CI.lower <- apply(h2Array, 1, quantile, (alpha / 2))
      
      ## CG EDITS (30 sept 19): Create data frame of communalities
      communalityDF <- data.frame(communalityDF, ## Obtained communalities
                                  h2SE,          ## Communality stand. errors
                                  h2CI.lower,    ## Lower CI bound
                                  h2CI.upper)    ## Upper CI bound
      
      ## CG EDITS (30 sept 19): Add column names to data frame of communalities
      colnames(communalityDF) <- 
        c("h2",  
          "SE", 
          paste0((alpha / 2) * 100, "th percentile"),
          paste0((1 - (alpha / 2)) * 100, "th percentile")) 
      
      #### ------- NAME BOOTSTRAP OBJECTS -------- ####
      
      ## Dimension names
      
      ## Name the factor indicators 
      rownames(loadSE) <-
        rownames(loadCI.upper) <-
        rownames(loadCI.lower) <- varNames
      
      ## Name the factors
      colnames(loadSE) <-
        colnames(loadCI.upper) <-
        colnames(loadCI.lower) <- 
        colnames(phiSE) <- rownames(phiSE) <-
        colnames(phiCI.upper) <- rownames(phiCI.upper) <-
        colnames(phiCI.lower) <- rownames(phiCI.lower) <-
        paste0("f", 1:numFactors)
      
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
        FICI.lower    <- 
        h2SE          <- 
        h2CI.upper    <- 
        h2CI.lower    <- NULL
    } # END if (bootstrapSE == FALSE)
    
    #### ------- BOOTSTRAP COMMUNALITIES -------- ####
    
    ## Check to ensure communalityDF and faH2 are equivalent
    CheckEquiv <- all.equal(communalityDF$h2, faXh2, 
                            check.attributes = FALSE, ## Ignoring attribute names
                            tolerance        = 1e-6)  ## precision for equivalence
    
    ## If unequivalent, give an ominous warning. 
    if (CheckEquiv == FALSE) {
      warning("Results may be untrustworthy. Unrotated communalities differ from rotated communalities.")
    } # END if (CheckEquiv == FALSE) 
    
    #### ------- NAME REMAINING OBJECTS -------- ####
    
    ## Designate the factor names
    facNames <- paste0("f", 1:numFactors)

    ## Name objects to have factor indicator names
    rownames(minLambda) <-
      ## CG EDITS (30 sept 19): Add row names to data frame of communalities
      rownames(communalityDF)   <- varNames
    
    ## Name objects to have factor names
    ## Change made March 4, 2019 NGW
    ## handle 1-fac models
    if(numFactors > 1) colnames(minLambda) <- facNames
    if(numFactors == 1) dimnames(minLambda) <- list(varNames, facNames)
    rownames(minPhi)  <-  
      colnames(minPhi)  <- facNames
    
    ## facIndeter is NA if urLoadings is specified
    if ( any(is.na(facIndeter)) == FALSE) rownames(facIndeter) <- facNames
    
    ## Add rotation convergence status to faFit
    faModelFit$convergedR <- uniqueSolutions[[1]]$RotationConverged
    
    #### ----- OUTPUT ----- ####
    
    fout <- list(R                     = R,
                 loadings              = minLambda,
                 Phi                   = minPhi,
                 facIndeterminacy      = facIndeter,
                 h2                    = communalityDF,
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
                 loadingsArray         = loadArray,
                 PhiArray              = phiArray,
                 facIndeterminacyArray = FIArray,
                 faControl             = faOut$faControl,
                 faFit                 = faModelFit,
                 rotate                = rotate,
                 rotateControl         = cnRotate,
                 unSpunSolution        = uniqueSolutions[[UnSpunPosition]],
                 targetMatrix          = targetMatrix,
                 Call                  = CALL)
    
    class(fout) <- 'faMain'
    fout
    
  } ## END faMAIN

