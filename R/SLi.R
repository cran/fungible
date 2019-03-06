#' Conduct a Schmid-Leiman Iterated (SLi) Target Rotation
#'
#' Compute an iterated Schmid-Leiman target rotation (SLi). This algorithm applies Browne's partially-specified Procrustes target rotation to obtain a full-rank bifactor solution from a rank-deficient (Direct) Schmid-Leiman procedure. Note that the target matrix is automatically generated based on the salient argument. Note also that the algorithm will converge when the partially-specified target pattern in the n-th iteration is equivalent to the partially-specified target pattern in the (n-1)th iteration.
#'
#' @param R (Matrix) A correlation matrix
#' @param SL (Matrix, NULL) A (rank-deficient) Schmid-Leiman (SL) bifactor solution (e.g., from a Schmid-Leiman or Direct Schmid-Leiman rotation). If NULL, the function will estimate the SL solution using the \code{\link{SchmidLeiman}} function.
#' @param rotate (Character) Designate which rotation algorithm to apply. See the \code{\link{faMain}} function for more details about possible rotations. A geomin rotation is the default.
#' @param numFactors (Vector) The number of latent factors at each level of analysis. For example, c(3, 1) estimates three latent factors in the first-order common factor model and one latent factor in the second-order common factor model (i.e., 3 group factors and 1 general factor).
#' @param salient (Numeric) A threshold parameter used to dichotomize factor loadings to create the target matrix. The default value is .20 (in absolute value) which is based on the Abad et al., 2017 application of this method.
#' @param urLoadings (Matrix, NULL) A full-rank matrix of unrotated factor loadings to be rotated using the (automatically generated) target matrix. If specified as NULL, a full-rank matrix of factor loadings will be extracted using the \code{\link{faX}} function. An unweighted least squares ("fals") procedure is the default.
#' @param freelyEstG (Logical) Specify whether the general factor loadings are freely estimated (in the partially-specified target matrix). If set to FALSE, only general factor loadings above the salient threshold will be estimated in the partially-specified target rotation.
#' @param gFac (Numeric, Vector) The position of the general factor(s) to be estimated. Solutions with multiple general factors may be estimated. Must either (a) freely estimate all loadings on the general factors or (b) only freely estimate general factor loadings that are above the salient threshold. The default column position is 1.
#' @param maxSLiItr (Numeric) The maximum number of iterations for the SLi procedure. Typically, 10 iterations is usually sufficient to converge (cf. Abad et al., 2017). The default is 20 iterations. 
#' @param digits (Numeric) The number of digits to round all output to. The default is no rounding.
#' @inheritParams faMain
#'
#' @return This function iterates the Schmid-Leiman target rotation and returns several relevant output.
#' \itemize{
#'   \item \strong{loadings}: (Matrix) The bifactor solution obtain from the SLi procedure.
#'   \item \strong{iterations}: (Numeric) The number of iterations required for convergence
#'   \item \strong{rotateControl}: (List) A list of the control parameters passed to the \code{\link{faMain}} function.
#'   \item \strong{faControl}: (List) A list of optional parameters passed to the factor extraction (\code{\link{faX}}) function.
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
#' @references Abad, F. J., Garcia-Garzon, E., Garrido, L. E., & Barrada, J. R. (2017). Iteration of partially specified target matrices: Application to the bi-factor case. \emph{Multivariate Behavioral Research, 52}(4), 416-429.
#' @references Giordano, C. & Waller, N. G. (under review). Recovering bifactor models: A comparison of seven methods.
#' @references Moore, T. M., Reise, S. P., Depaoli, S., & Haviland, M. G. (2015). Iteration of partially specified target matrices: Applications in exploratory and Bayesian confirmatory factor analysis. \emph{Multivariate Behavioral Research, 50}(2), 149-161.
#' @references Reise, S. P., Moore, T. M., & Haviland, M. G. (2010). Bifactor models and rotations: Exploring the extent to which multidimensional data yield univocal scale scores. \emph{Journal of Personality Assessment, 92}(6), 544-559.
#' @references Schmid, J., & Leiman, J. M. (1957). The development of hierarchical factor solutions. \emph{Psychometrika, 22}(1), 53-61.
#'
#' @examples
#' ## Generate a bifactor model
#' bifactor <- matrix(c(.35, .61, .00, .00,
#'                      .35, .61, .00, .00,
#'                      .35, .61, .00, .00,
#'                      .35, .00, .61, .00,
#'                      .35, .00, .61, .00,
#'                      .35, .00, .61, .00,
#'                      .35, .00, .00, .61,
#'                      .35, .00, .00, .61,
#'                      .35, .00, .00, .61),
#'                    nrow = 9, ncol = 4, byrow = TRUE)
#'
#' ## Model-implied correlation (covariance) matrix
#' R <- bifactor %*% t(bifactor)
#'
#' ## Unit diagonal elements
#' diag(R) <- 1
#'
#' Out1 <- SLi(R          = R,
#'             numFactors = c(3, 1),
#'             digits     = 2)
#'
#' @export
#'

SLi <-
  function(R,                         ## Correlation matrix
           SL            = NULL,      ## SL loading matrix
           rotate        = "geominQ", ## Oblique rotation for computing SL
           numFactors    = NULL,      ## Number of SL group factors to estimate
           facMethod     = "fals",     ## Factor extraction method
           salient       = .20,       ## Salient threshold
           urLoadings      = NULL,      ## Unrotated factor structure
           freelyEstG    = TRUE,      ## Freely estimate all general factor loadings?
           gFac          = 1,         ## Column(s) of general factor(s)
           maxSLiItr     = 20,        ## Max num of iterations
           digits        = NULL,      ## Round the output
           rotateControl = NULL,      ## Control rotation tuning parameters
           faControl     = NULL){     ## Control parameters of factor extraction
    
    ## ~~~~~~~~~~~~~~~~~~ ##
    #### Error Checking ####
    ## ~~~~~~~~~~~~~~~~~~ ##
    
    ## Is correlation matrix symmetric
    if ( !isSymmetric(R) ) {
      stop("The user-defined correlation matrix is not symmetric")
    } # END if (!isSymmetric(R)) {
    
    ## IF SL is to be estimated, give the correct argument inputs.
    if ( is.null(SL) ) {
      
      if ( is.null(numFactors) ) {
        stop("To estimate the SL solution, please specify the number of factors at each level of the higher-order solution.")
      } # END if ( is.null(numFactors) )
      
      ## Number of factors correctly specified?
      if (length(numFactors) <= 1) {
        stop("The 'numFactors' argument must be a vector of 2 or 3 values.")
      } # END if (length(numFactors) <=1)
      
    } # END if ( is.null(SL) )
    
    
    ## Must specify either R or urLoadings
    if ( is.null(urLoadings) & is.null(R) ) {
      stop("Must specify either a correlation matrix (for factor extraction) or an unrotated factor structure.")
    } # END if ( is.null(urLoadings) & is.null(R) )
    
    # ## Check rotateControl
    # 
    # ## rotateControl specified?
    # if ( is.null(rotateControl) ) {
    # 
    #   ## Set the default rotateControl values (before user-specification of values)
    #   cnDefault <- list(gamma       = 0,
    #                     delta       = .01,
    #                     kappa       = 0,
    #                     k           = NULL,
    #                     standardize = "none",
    #                     epsilon     = 1e-5,
    #                     power       = 4,
    #                     maxItr      = 15000)
    # 
    #   ## Used as func output, if it is NULL, need to specify
    #   ## Full error checking takes place in faMain()
    #   rotateControl <- cnDefault
    # } # END if ( is.null(rotateControl) )
    
    # ## Check faControl
    # 
    # ## If faControl is not specified, give it the defaults (used for func output)
    # if ( is.null(faControl) ) {
    # 
    #   ## Set the default values of all control agruments
    #   cnFA <- list(treatHeywood   = TRUE,
    #                nStart         = 10,
    #                maxCommunality = .995,
    #                epsilon        = 1e-4,
    #                communality    = "SMC",
    #                maxItr         = 15000)
    # 
    #   ## Used as func output, if it is NULL, need to specify
    #   ## Full error checking takes place in faX()
    #   faControl <- cnFA
    # 
    # } # END if ( is.null(faControl) )
    
    ## If digits is not specified, give it an arbitrarily large value
    if ( is.null(digits) ) {
      digits <- options()$digits
    } # END if ( is.null(digits) )
    
    ## ~~~~~~~~~~~~~~ ##
    #### Begin Code ####
    ## ~~~~~~~~~~~~~~ ##
    
    ## If the SL solution is not provided, estimate it
    if ( is.null(SL) ) {
      
      ## Call the SL function and return the rank-deficient factor solution
      SL <- SchmidLeiman(R             = R,
                         numFactors    = numFactors,
                         rotate        = rotate,
                         facMethod     = facMethod,
                         rotateControl = rotateControl,
                         faControl     = faControl)$B
      
    } # END if ( is.null(SL) )
    
    ## Number of TOTAL factors
    numberFactors <- ncol(SL)
    
    ## Generate a target matrix by specifying values below the salient
    ## threshold to zero and allowing all NAs to be freely estimated
    targetOld <- SL
    targetOld[abs(targetOld) >= salient] <- NA
    targetOld[abs(targetOld) <  salient] <- 0
    
    ## Is the general factor freely estimated?
    if (freelyEstG == TRUE) {
      
      ## Set all general factor loadings to NA in order to be freely estimated
      targetOld[, gFac] <- NA
      
    } # END if (freelyEstG == TRUE)
    
    ## Conduct new factor analysis on numberFactors factors
    if ( is.null(urLoadings) ) {
      
      ## Extract a full-rank factor structure matrix
      faOut <-
        faX(R          = R,
            numFactors = numberFactors,
            facMethod  = facMethod,
            faControl  = faControl)
      
      ## Save the unrotated factor loadings
      fact.loadings <- faOut$loadings[]
      
      ## Check for Heywood case
      if ( faOut$faFit$Heywood == TRUE ) {
        warning("There is a Heywood case in the unrotated, full-rank factor structure matrix.")
      } # END if ( faOut$faFit$Heywood == TRUE )
      
    } else {
      
      ## Use the pre-specified unrotated factor structure
      fact.loadings <- urLoadings
      
    } # END if (is.null(urLoadings)) {
    
    ## Begin a count for the number of iterations used for convergence
    numIters <- 0
    
    ## Iterate until convergence OR hitting maxSLiItr limit
    for (iter in 1:maxSLiItr){
      
      ## Procrustes rotation using the targetOld patterned matrix
      SLt <-   faMain(urLoadings    = fact.loadings,
                      rotate        = "targetT",
                      targetMatrix  = targetOld, ## This is updated each loop
                      rotateControl = rotateControl)
      
      ## Generate new specified target matrix
      targetNew <- SLt$loadings[]
      
      ## Set loadings to be freely estimated
      targetNew[abs(targetNew) >= salient] <- NA
      
      ## Set loadings set to zero
      targetNew[abs(targetNew) <  salient] <- 0
      
      ## Is the general factor freely estimated or only salient loadings?
      if (freelyEstG == TRUE) {
        
        ## Set all general factor loadings to be freely estimated
        targetNew[, gFac] <- NA
        
      } # END if (freelyEstG == TRUE)
      
      ## If old and new target matrix specification are perfectly the same, stop
      if (all.equal(targetNew, targetOld, check.attributes = FALSE) == TRUE) {
        
        ## break will stop the for loop
        break
        
      } else {
        
        ## Rename for the next loop in the iteration
        targetOld <- targetNew
        
        ## Use previously rotated loadings as the new matrix to preform Procrustes on
        fact.loadings <- SLt$loadings
        
        ## Add +1 to number of iterations
        numIters <- numIters + 1
        
      } # END if (all.equal(targetNew, targetOld, check.attributes = FALSE) == TRUE)
      
    } # End the for loop of iterations
    
    ## Display note if the iteration did not converge
    if (iter == maxSLiItr & iter > 1) {
      
      warning("The solution did not converge in the maximum allotted iterations.
              \n Try increasing the maxSLiItr argument.")
      
    } # END if (iter == maxSLiItr & iter > 1) {
    
    ## Extract the final estimated solution
    SLi.loadings <- SLt$loadings
    
    ## Give proper column names
    colnames(SLi.loadings) <-
      c(paste(rep("g", numFactors[2]), c(1:numFactors[2]), sep = ""),
        paste(rep("F", numFactors[1]), c(1:numFactors[1]), sep = ""))
    
    ## Reflect all factors by default
    SLiSign <- diag( sign( colSums(SLi.loadings) ) )
    SLi.loadings <- SLi.loadings %*% SLiSign
    
    ## Final output list
    ##    [[1]]: the final iterated Schmid-Leiman matrix of loadings
    ##    [[2]]: number of iterations to achieve final output
    
    list(loadings      = round(SLi.loadings, digits),
         iterations    = round(numIters, digits),
         rotateControl = rotateControl,
         faControl     = faControl)
    
  } # END SchmidLeimanIterated()
