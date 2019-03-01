#' Automatic Factor Rotation from Random Configurations with Bootstrap Standard Errors
#' 
#' This function conducts factor rotations (using the \pkg{GPArotation} package) 
#' from a user-specified number of random (orthogonal) starting configurations. 
#' Based on the resulting discrepancy function, the function determines the 
#' number of local minima and, among these local solutions, will find the 
#' "global minimum" (i.e., the minimized discrepancy value from the finite 
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
#' factoring, and "pca" for principal components analysis. The default method 
#' is "fals". 
#' \itemize{
#'   \item \strong{"fals"}: Factors are extracted using the unweighted least 
#'   squares estimation procedure using the \code{\link{fals}} function.
#'   \item \strong{"faml"}: Factors are extracted using the maximum likelihood 
#'   estimation procedure using the \code{\link[stats]{factanal}} function.
#'   \item \strong{"fapa"}: Factors are extracted using the iterated principal 
#'   axis factoring estimation procedure using the \code{\link{fapa}} function.
#'   \item \strong{"pca}: Principal components are extracted. 
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
#' @param targetMatrix (Matrix) The target matrix for (fully and partially) 
#' specified target rotations. To conduct Browne's (2001) partially-specified 
#' target rotation, freely estimated factor loadings are designated by "NA" values.
#' @param bootstrapSE (Logical) Computes bootstrap standard errors. All bootstrap 
#' samples are aligned to the global minimum solution. Defaults to 
#' bootstrapSE = FALSE (no standard errors). 
#' @param numBoot (Numeric) The number bootstraps. Defaults to numBoot = 1000.
#' @param CILevel (Numeric) The confidence level (between 0 and 1) of the bootstrap 
#' confidence interval. Defaults to CILevel = .95.
#' @param Seed (Numeric) Starting seed for reproducible bootstrap results. 
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
#'     \item \strong{"maxRsqr"}: Initial communalities equal the largest squared 
#'     correlation in each column of the correlation matrix.
#'     \item \strong{"unity"}: Initial communalities equal 1.0 for all variables.
#'   }
#'   \item \strong{maxItr}: (Numeric) In \code{fapa}, the maximum number of 
#'   iterations to reach convergence. Defaults to maxItr = 15,000.
#' }
#' 
#' @param rotateControl (List) A list of control values to pass to the factor rotation algorithms.
#' \itemize{
#'   \item \strong{numberStarts}: (Numeric) The number of random (orthogonal) 
#'   starting configurations for the chosen rotation method (e.g., oblimin). 
#'   Defaults to numberStarts = 10. 
#'   \item \strong{itemSort}: (Logical) If TRUE, sort the row order of all the following output 
#'   such that variables loading on a common factor are grouped together for 
#'   ease of interpretation: (a) the global minimum factor loadings, 
#'   (b) indicator communalities, (c) factor-loading bootstrap standard errors, 
#'   (d) factor-loading bootstrap confidence interval quantiles (both upper and 
#'   lower), and (e) the array of all factor-loading bootstrap results. 
#'   Defaults to itemSort = FALSE.
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
#'   rotation discrepancy values (to the number of  digits specified in the \code{epsilon} 
#'   argument of the \code{rotateControl} list) into sets with equivalent values. For example, 
#'   by default \code{epsilon = 1e-5.} and thus \code{} will only evaluate the discrepancy 
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
#'   in ascending order of their discrepancy values (i.e., the first solution 
#'   is the "global" minimum). Each solution returns the 
#'   (a) factor loadings, 
#'   (b) factor correlations, 
#'   (c) the discrepancy value of the rotation algorithm, 
#'   (d) A vector of factor indeterminacy indices for each common factor, and 
#'   (d) whether the rotation procedure converged.
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
#'     \item \strong{converged}: (Logical) TRUE if the least squares or 
#'     principal axis factor extraction routine converged. 
#'   }
#'   \item \strong{rotateControl}: (List) A list of the control parameters 
#'   passed to the rotation algorithm.
#'   \item \strong{itemOrder}: (Vector) The final item order if \code{itemSort = TRUE}. 
#'   \item \strong{Call}: (call) A copy of the function call.
#' }
#' 
#' @references Browne, M. W. (1972). Orthogonal rotation to a partially specifed
#' target. \emph{British Journal of Statistical Psychology, 25}(1), 115-120.
#' @references Browne, M. W. (2001). An overview of analytic rotation in 
#' exploratory factor analysis. \emph{Multivariate Behavioral Research, 36}(1), 111-150.
#' @references Cureton, E. E., & Mulaik, S. A. (1975). The weighted varimax 
#' rotation and the promax rotation. \emph{Psychometrika, 40}(2), 183-195.
#' @references Guttman, L. (1955). The determinacy of factor score matrices with 
#' implications for five other basic problems of common factor theory. 
#' \emph{British Journal of Statistical Psychology, 8}(2), 65-81.
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
#'   \item Casey Giordano (Giord023@umn.edu)
#'   \item Niels G. Waller (nwaller@umn.edu)
#'   \item The authors thank Allie Cooperman and Hoang
#'    Nguyen for their help  implementing the standard error estimation and the 
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
#'                                     itemSort     = TRUE,
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
           targetMatrix  = NULL,      ## If targetT, matrix of specified loadings
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
    
    ## Check rotation
    
    ## That a rotation is specified
    if ( !is.character(rotate) ) {
      stop("The 'rotate' argument must be a character string. See ?rotate for available options.")
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
      
      ## Class must be matrix or DR to analyze via cor() function
      if ( class(X) %in% c("matrix", "data.frame", "loadings") == FALSE) {
        stop("'X' must be of class matrix, data.frame, or loadings.")
      } # END if ( class(X) %in% c("matrix", "data.frame", "loadings") == FALSE)
      
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
                     itemSort     = FALSE,
                     gamma        = 0,
                     delta        = .01,
                     kappa        = 0,
                     k            = nrow(urLoadings),
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
    internalRotate <- function(lambda, rotation, rotateControl) {
      ## Purpose: Simple rotation wrapper
      ##
      ## Args: lambda:        (matrix) unrotated loadings to rotate
      ##       rotation:      (character) which rotation to use
      ##       rotateControl: (list) tuning parameters to pass to rotation
      ##
      
      ## Determine column dimensions
      matrixDim <- ncol(lambda)
      
      ## Perform rotation with the specified parameters
      switch(rotation,
             "none" = {
               list(loadings = lambda,
                    Phi      = NULL)
             },
             "oblimin" = {
               GPArotation::oblimin(lambda,
                                    Tmat      = randStart(matrixDim),
                                    gam       = cnRotate$gamma,
                                    normalize = FALSE,
                                    eps       = cnRotate$epsilon,
                                    maxit     = cnRotate$maxItr)
             },
             "quartimin" = {
               GPArotation::quartimin(lambda,
                                      Tmat      = randStart(matrixDim),
                                      maxit     = cnRotate$maxItr,
                                      eps       = cnRotate$epsilon,
                                      normalize = FALSE)
             },
             "targetT" = {
               GPArotation::targetT(lambda,
                                    Tmat      = randStart(matrixDim),
                                    Target    = targetMatrix,
                                    normalize = FALSE,
                                    eps       = cnRotate$epsilon,
                                    maxit     = cnRotate$maxItr)
             },
             "targetQ" = {
               GPArotation::targetQ(lambda,
                                    Tmat      = randStart(matrixDim),
                                    Target    = targetMatrix,
                                    normalize = FALSE,
                                    eps       = cnRotate$epsilon,
                                    maxit     = cnRotate$maxItr)
             },
             "oblimax" = {
               GPArotation::oblimax(lambda,
                                    Tmat      = randStart(matrixDim),
                                    normalize = FALSE,
                                    eps       = cnRotate$epsilon,
                                    maxit     = cnRotate$maxItr)
             },
             "entropy" = {
               GPArotation::entropy(lambda,
                                    Tmat      = randStart(matrixDim),
                                    normalize = FALSE,
                                    eps       = cnRotate$epsilon,
                                    maxit     = cnRotate$maxItr)
             },
             "quartimax" = {
               GPArotation::quartimax(lambda,
                                      Tmat      = randStart(matrixDim),
                                      normalize = FALSE,
                                      eps       = cnRotate$epsilon,
                                      maxit     = cnRotate$maxItr)
             },
             "varimax" = {
               GPArotation::Varimax(lambda,
                                    Tmat      = randStart(matrixDim),
                                    normalize = FALSE,
                                    eps       = cnRotate$epsilon,
                                    maxit     = cnRotate$maxItr)
             },
             "simplimax" = {
               GPArotation::simplimax(lambda,
                                      Tmat      = randStart(matrixDim),
                                      normalize = FALSE,
                                      eps       = cnRotate$epsilon,
                                      k         = cnRotate$k,
                                      maxit     = cnRotate$maxItr)
             },
             "bentlerT" = {
               GPArotation::bentlerT(lambda,
                                     Tmat      = randStart(matrixDim),
                                     maxit     = cnRotate$maxItr,
                                     eps       = cnRotate$epsilon,
                                     normalize = FALSE)
             },
             "bentlerQ" = {
               GPArotation::bentlerQ(lambda,
                                     Tmat      = randStart(matrixDim),
                                     maxit     = cnRotate$maxItr,
                                     eps       = cnRotate$epsilon,
                                     normalize = FALSE)
             },
             "tandemI" = {
               GPArotation::tandemI(lambda,
                                    Tmat      = randStart(matrixDim),
                                    maxit     = cnRotate$maxItr,
                                    eps       = cnRotate$epsilon,
                                    normalize = FALSE)
             },
             "tandemII" = {
               GPArotation::tandemII(lambda,
                                     Tmat      = randStart(matrixDim),
                                     maxit     = cnRotate$maxItr,
                                     eps       = cnRotate$epsilon,
                                     normalize = FALSE)
             },
             "geominT" = {
               GPArotation::geominT(lambda,
                                    Tmat      = randStart(matrixDim),
                                    delta     = cnRotate$delta,
                                    normalize = FALSE,
                                    eps       = cnRotate$epsilon,
                                    maxit     = cnRotate$maxItr)
             },
             "geominQ" = {
               GPArotation::geominQ(lambda,
                                    Tmat      = randStart(matrixDim),
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
                                Tmat      = randStart(matrixDim),
                                kappa     = cnRotate$kappa,
                                maxit     = cnRotate$maxItr,
                                eps       = cnRotate$epsilon,
                                normalize = FALSE)
             },
             "cfQ" = {
               GPArotation::cfQ(lambda,
                                Tmat      = randStart(matrixDim),
                                kappa     = cnRotate$kappa,
                                eps       = cnRotate$epsilon,
                                normalize = FALSE,
                                maxit     = cnRotate$maxItr)
             },
             "infomaxT" = {
               GPArotation::infomaxT(lambda,
                                     Tmat      = randStart(matrixDim),
                                     normalize = FALSE,
                                     eps       = cnRotate$epsilon,
                                     maxit     = cnRotate$maxItr)
             },
             "infomaxQ" = {
               GPArotation::infomaxQ(lambda,
                                     Tmat      = randStart(matrixDim),
                                     normalize = FALSE,
                                     eps       = cnRotate$epsilon,
                                     maxit     = cnRotate$maxItr)
             },
             "mccammon" = {
               GPArotation::mccammon(lambda,
                                     Tmat      = randStart(matrixDim),
                                     normalize = FALSE,
                                     eps       = cnRotate$epsilon,
                                     maxit     = cnRotate$maxItr)
             },
             "bifactorT" = {
               GPArotation::bifactorT(lambda,
                                      Tmat      = randStart(matrixDim),
                                      normalize = FALSE,
                                      eps       = cnRotate$epsilon,
                                      maxit     = cnRotate$maxItr)
             },
             "bifactorQ" = {
               GPArotation::bifactorQ(lambda,
                                      Tmat      = randStart(matrixDim),
                                      normalize = FALSE,
                                      eps       = cnRotate$epsilon,
                                      maxit     = cnRotate$maxItr)
             })
    } # END internalRotate
    
    ## Function to re-order the factor loading/correlation matrices
    orderFactors <- function(Lambda,
                             PhiMat) {
      ## Purpose: Sort the factors into conventional order
      ##
      ## Args: Lambda: (Matrix) Factor loading matrix
      ##       PhiMat: (Matrix) Factor correlation matrix

      ## Wipe out low loadings (semi-arbitrarily)
      F0 <- Lambda
      F0[abs(Lambda) < .3] <- 0

      ## Flip poorly keyed factors (and associated factor correlations)
      FsgnMat <- diag(sign(apply(F0, 2, sum)))
      Lambda  <- Lambda %*% FsgnMat
      Phi     <- FsgnMat %*% PhiMat %*% FsgnMat

      ## Vector used to re-order the factor loading/correlation matrices
      ## Based on sum of squared loadings
      Lambda.SS <- sort.list( apply(F0^2, 2, sum), decreasing = TRUE )

      ## Sort the keyed factor loadings
      Lambda <- Lambda[, Lambda.SS]

      ## Sort the keyed factor correlations
      I   <- diag(ncol(Lambda))
      Phi <- I[Lambda.SS, ] %*% Phi %*% I[, Lambda.SS]

      ## Return a list of the output
      list(Lambda = Lambda,
           PhiMat = Phi)

    }  #END orderFactors
    
    ## Compute Guttman's factor determinacy indices
    GuttmanIndices <- function(Lambda, 
                               PhiMat) {
      ## Purpose: Compute Guttman (1955) factor indeterminacy indices
      ##
      ## Args:    Lambda: (Matrix) Rotated factor loadings matrix
      ##          PhiMat: (Matrix) Factor correlation matrix
      
      ## Fator structure (works for either oblique or orthogonal model)
      facStruct <- Lambda %*% PhiMat
      
      ## Model-implied correlation matrix
      modImpR <- facStruct %*% t(Lambda)
      
      ## Otherwise, diagonal elements are the communalities
      diag(modImpR) <- 1
      
      ## Factor indeterminacy solution
      sqrt( diag( t(facStruct) %*% solve(modImpR) %*% facStruct))
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
      
      ## No model fit stuff
      faModelFit <- NULL
      
      ## Compute communalities before rotation
      faXh2 <- apply(urLoadings^2, 1, sum)
      
    } # END if ( !is.null(urLoadings) ) 
    
    ## If urLoadings is not supplied, extract it
    if ( is.null(urLoadings) ) {
      
      ## Call faX function to extract factor structure matrix
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
      faModelFit <- faOut$faFit
      
    } # END if ( is.null(urLoadings) )
    
    ## If urLoadings is provided and numFactors isn't, specify numFactors
    if ( is.null(numFactors) && !is.null(urLoadings) ) {
      numFactors <- ncol(urLoadings)
    } # END if ( is.null(numFactors) && !is.null(urLoadings) )
    
    #### ------- STANDARDIZATION -------- ####
    
    ## Call standardize function
    Stnd <- Standardize(method = cnRotate$standardize,
                        lambda = urLoadings)
    
    ## Extract DvInv for later unstandardization
    lambda <- Stnd$lambda
    
    #### ------- ROTATION -------- ####
    
    ## Pre-allocate a list for the different attempts
    starts <- vector("list",
                     cnRotate$numberStarts)
    
    ## For each list element, do the specified rotate and save all output
    starts <- lapply(starts, function(x) {
      
      ## Rotate the standardized factor structure matrix
      rotatedLambda <- internalRotate(lambda        = lambda,
                                      rotation      = rotate,
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
    
    #### ------- COUNT UNIQUE SOLUTIONS -------- ####
    
    ## Save the minimized discrepancy function for each attempt
    ##    $Table[,2] is the value of criterion at each iteration, which.min
    ##    finds the minimum criterion value across each attempt
    
    ## Evaluate the minimum disc functions to find smallest value
    if (rotate == "promaxQ") {
      
      ## For each random start, find the evaluated discrepancy function
      DiscrepFunc <- sapply(starts, function(attempt) attempt$vmaxDiscrepancy)
      
    } else {
      
      ## For each random start, find the evaluated discrepancy function
      DiscrepFunc <- sapply(starts, function(attempt) min(attempt$Table[, 2]))
      
    } # END if (rotate == "promaxQ")
    
    ## Find the sort order (minimum first)
    sortOrder <- order(DiscrepFunc)
    
    ## Create new list to hold the output
    uniqueSolutions <- vector("list", length(sortOrder))
    
    ## Create list of output from local solutions
    for (numStart in 1:length(sortOrder)) {
      
      ## Determine which element of sortOrder to grab
      num <- which(sortOrder == numStart)
      
      ## Extract the relevant factor loadings and Phi matrices
      selected <- starts[[num]]
      
      #### ------- FACTOR SORT -------- ####
      
      ## sort loadings and Phi matrix
      sortedSols <- orderFactors(Lambda = selected$loadings,
                                 PhiMat = selected$Phi)
      
      ## Overwrite "selected" list
      selected$loadings <- sortedSols$Lambda
      selected$Phi      <- sortedSols$PhiMat
      
      ## Select the factor loadings
      uniqueSolutions[[numStart]]$loadings <- selected$loadings
      
      ## Select factor correlations
      uniqueSolutions[[numStart]]$Phi <- selected$Phi
      
      ## Select discrepancy values
      uniqueSolutions[[numStart]]$DiscrepValue <- DiscrepFunc[num]
      
      #### ----- FACTOR INDETERMINACY ----- ####
      
      ## Guttman's factor indeterminacy indices
      uniqueSolutions[[numStart]]$facIndeterminacy <- 
        GuttmanIndices(Lambda = selected$loadings, 
                       PhiMat = selected$Phi)
      
      ## Did the local optima solution converge
      uniqueSolutions[[numStart]]$converged <- selected$convergence
      
    } # END for (numStart in 1:length(sortOrder))
    
    ## Extracted discrepancy values from uniqueSolutions to find solution sets
    DisVal <- unlist(lapply(uniqueSolutions, function(x) x$DiscrepValue))
    
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
    facIndeter <- uniqueSolutions[[1]]$facIndeterminacy
    
    #### -------- BOOTSTRAP SETUP -------- ####
    
    ## If true, compute bootstrap standard errors
    if (bootstrapSE == TRUE) {
      
      # Number of scale items
      nVar <- ncol(X)
      
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
      
      #### -------- BOOTSTRAP FOR LOOP -------- ####
      
      ## Analyses on 'numBoot' number of random resamples of X matrix
      for (iSample in seq_len(numBoot)) {
        
        ## Set the seed for reproducibility
        set.seed(iSample + Seed)
        
        ## Resample (with replacement) from X raw data matrix
        bsSample <- sample(x       = rows,
                           size    = nSubj,
                           replace = TRUE)
        
        ## Create correlation matrix
        Rsamp <- cor(X[bsSample, ], ...)
        
        ## Extract unrotated factors using resampled data matrix
        bsLambda <- faX(R          = Rsamp,
                        numFactors = numFactors,
                        facMethod  = facMethod,
                        faControl  = faControl)$loadings
        
        ## Conduct standardization
        bsStnd <- Standardize(method = cnRotate$standardize,
                              lambda = bsLambda)
        
        ## Extract the standardized bootstrapped (unrotated) factor loadings
        bsLambda <- bsStnd$lambda
        
        ## Find the "global" min out of all the random starting configurations
        bsStarts <- vector("list", cnRotate$numberStarts)
        
        ## Conduct rotations from random start values
        bsStarts <- lapply(bsStarts, function(x) {
          
          ## Rotate the bootstrapped samples
          bsRotated <- internalRotate(lambda        = bsLambda,
                                      rotation      = rotate,
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
        if (rotate == "promaxQ") {
          
          ## For each random start, find the evaluated discrepancy function
          bsDiscrepFunc <- sapply(bsStarts, function(attempt) attempt$vmaxDiscrepancy)
          
        } else {
          
          ## For each random start, find the evaluated discrepancy function
          bsDiscrepFunc <- sapply(bsStarts, function(attempt) min(attempt$Table[, 2]))
          
        } # END if (rotate == "promaxQ")
        
        ## Of all random configs, determinine which has the lowest criterion value
        bsLambda <- bsStarts[[which.min(DiscrepFunc)]]$loadings
        bsPhi    <- bsStarts[[which.min(DiscrepFunc)]]$Phi
        
        ## Unstandardize the minimum rotated solution for this bootstrap sample
        bsLambda <- bsStnd$DvInv %*% bsLambda
        
        ## Align bootstrap sample with global minimum from before
        Aligned <- faAlign(F1   = minLambda,
                           F2   = bsLambda,
                           Phi2 = bsPhi)
        
        ## Save loadings as the bootstrapped factor loadings
        fList[[iSample]]   <- Aligned$F2
        
        ## Save correlations as the bootstrapped factor correlations
        phiList[[iSample]] <- Aligned$Phi2
        
        ## Save factor indeterminacy indices
        FIList[[iSample]] <- GuttmanIndices(Lambda = Aligned$F2,
                                            PhiMat = Aligned$Phi2)
        
      } # END for (iSample in seq_len(numBoot))
      
      #### ----- STANDARD ERRORS ----- ####
      
      # Convert list of matrices into array
      loadArray <- array(unlist(fList), c(nVar, numFactors, numBoot))
      phiArray  <- array(unlist(phiList), c(numFactors, numFactors, numBoot))
      FIArray   <- array(unlist(FIList), c(1, numFactors, numBoot))
      
      
      # Bootstrap standard errors for factor loadings
      loadSE <- apply(loadArray, 1:2, sd)
      loadSE <- round(loadSE, digits)
      
      # Bootstrap standard errors for factor correlations
      phiSE <- apply(phiArray, 1:2, sd)
      phiSE <- round(phiSE, digits)
      
      ## Bootstrap standard errors for factor indeterminacy
      FISE <- apply(FIArray, 1:2, sd)
      FISE <- round(FISE, digits)
      
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
        colnames(FISE) <- paste0("f", 1:numFactors)
      
    } # END if (bootstrapSE == TRUE)
    
    ## If not specified, return NULL
    if (bootstrapSE == FALSE) {
      loadSE         <-
        loadCI.upper <-
        loadCI.lower <-
        loadArray    <-
        phiSE        <-
        phiCI.upper  <-
        phiCI.lower  <-
        phiArray     <- 
        FIArray      <- 
        FISE         <- NULL
    } # END if (bootstrapSE == FALSE)
    
    #### ------- COMMUNALITIES -------- ####
    
    ## Compute communalities. If orthogonal model, Phi is identity matrix
    rotateH2 <- diag(minLambda %*% minPhi %*% t(minLambda))
    
    ## Check to ensure rotateH2 and faH2 are equivalent
    CheckEquiv <- all.equal(rotateH2, faXh2, 
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
      names(rotateH2)   <- varNames
    
    ## Name objects to have factor names
    colnames(minLambda) <- 
      rownames(minPhi)  <-  
      colnames(minPhi)  <- 
      names(facIndeter) <- facNames
    
    #### ------- SORT ITEMS -------- ####
    
    if (cnRotate$itemSort == TRUE) {
      
      ## Save all output from item sorting function
      itemSorting <- faSort(fmat = minLambda,
                            phi  = minPhi)
      
      ## Determine the order in which items are sorted
      sortOrder <- itemSorting$sortOrder
      
      ## Sort items in global min lambda
      minLambda <- minLambda[sortOrder, ]

      ## Sort items in rotateH2 (communality of indicators)
      rotateH2 <- rotateH2[sortOrder]

      ## If bootstrap is done, re-order factor indicators
      if (bootstrapSE == TRUE) {
        loadSE       <- loadSE[sortOrder, ]
        loadCI.upper <- loadCI.upper[sortOrder, ]
        loadCI.lower <- loadCI.lower[sortOrder, ]
        loadArray    <- loadArray[sortOrder, , ]
        
      } # END if (bootstrapSE == TRUE)
      
    } # END if (cnRotate$itemSort == TRUE) 
    
    #### ----- OUTPUT ----- ####
    
    list(R                     = R,
         loadings              = round(minLambda, digits),
         Phi                   = round(minPhi,    digits),
         facIndeterminacy      = facIndeter,
         h2                    = rotateH2,
         loadingsSE            = loadSE,
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
         faControl             = faControl,
         faFit                 = faModelFit,
         rotateControl         = cnRotate,
         itemOrder             = sortOrder,
         Call                  = CALL)
    
  } ## END rotate()
