#' Conduct a Schmid-Leiman Iterated (SLi) Target Rotation
#'
#' Compute an iterated Schmid-Leiman target rotation (SLi). This algorithm applies Browne's partially-specified Procrustes target rotation to obtain a full-rank bifactor solution from a rank-deficient (Direct) Schmid-Leiman procedure. Note that the target matrix is automatically generated based on the salient argument. Note also that the algorithm will converge when the partially-specified target pattern in the n-th iteration is equivalent to the partially-specified target pattern in the (n-1)th iteration.
#'
#' @param R (Matrix) A correlation matrix
#' @param SL (Matrix, NULL) A (rank-deficient) Schmid-Leiman (SL) bifactor solution (e.g., from a Schmid-Leiman or Direct Schmid-Leiman rotation). If NULL, the function will estimate the SL solution using the \code{\link{SchmidLeiman}} function.
#' @param rotate (Character) Designate which rotation algorithm to apply. See the \code{\link{rotate}} function for more details about possible rotations. A geomin rotation is the default.
#' @param numFactors (Vector) The number of latent factors at each level of analysis. For example, c(3, 1) estimates three latent factors in the first-order common factor model and one latent factor in the second-order common factor model (i.e., 3 group factors and 1 general factor).
#' @param facMethod (Character) The method used for factor extraction (see \code{\link{faX}} for more details). The supported options are "fals" for unweighted least squares, "faml" for maximum likelihood, and "fapa" for iterated principal axis factoring. The default method is "fals".
#' \itemize{
#'   \item \strong{"fals"}: Factors are extracted using the unweighted least squares estimation procedure using the \code{\link[fungible]{fals}} function.
#'   \item \strong{"faml"}: Factors are extracted using the maximum likelihood estimation procedure using the \code{\link[stats]{factanal}} function.
#'   \item \strong{"fapa"}: Factors are extracted using the iterated principal axis factoring estimation procedure using the \code{\link{fapa}} function.
#' }
#' @param salient (Numeric) A threshold parameter used to dichotomize factor loadings to create the target matrix. The default value is .20 (in absolute value) which is based on the Abad et al., 2017 application of this method.
#' @param urLambda (Matrix, NULL) A full-rank matrix of unrotated factor loadings to be rotated using the (automatically generated) target matrix. If specified as NULL, a full-rank matrix of factor loadings will be extracted using the \code{\link{faX}} function. An unweighted least squares ("fals") procedure is the default.
#' @param freelyEstG (Logical) Specify whether the general factor loadings are freely estimated (in the partially-specified target matrix). If set to FALSE, only general factor loadings above the salient threshold will be estimated in the partially-specified target rotation.
#' @param gFac (Numeric, Vector) The position of the general factor(s) to be estimated. Solutions with multiple general factors may be estimated. Must either (a) freely estimate all loadings on the general factors or (b) only freely estimate general factor loadings that are above the salient threshold. The default column position is 1.
#' @param numberStarts (Numeric) The number of random starting configurations used for each iteration of the partially-specified target rotation. Note that if the algorithm requires 5 iterations for convergence, 10 random starts will require a total of 50 partially-specified target rotations. The default value is 10.
#' @param maxITR (Numeric) The maximum number of iterations attempted before failing to converge. Typically, 10 iterations is usually sufficient to converge (cf. Abad et al., 2017).
#' @param digits (Numeric) The number of digits to round all output to. The default is no rounding.
#' @param rotateControl (List) A list of control values to pass to the factor rotations in the \code{\link{rotate}} function.
#' \itemize{
#'   \item \strong{gamma}: (Numberic) This is a tuning parameter (between 0 and 1, inclusive) for an oblimin rotation.  See the GPArotation library's oblimin documentation for more details.
#'   \item \strong{delta}: (Numberic) This is a tuning parameter for the geomin rotation. It adds a small number (default = .01) to factor loadings before computing the geometric means in the discrepancy function.
#'   \item \strong{kappa}: (Numeric) The main parameterization of the Crawford-Ferguson family of factor rotations. Specific values correspond to certain rotations. For instance (adapted from GPArotation's help page), 0=Quartimax, 1/p=varimax, m/(2*p)=Equamax, (m-1)/(p+m-2)=Parsimax, 1=Factor parsimony. Note that 'p' is the number of factor indicators and 'm' is the number of common factors. The default value is 0 (Quartimax).
#'   \item \strong{k}: (Numeric) A specific parameter of the simplimax rotation. The default value is nrow(lambda).
#'   \item \strong{epsilon}: (Numeric) The rotational convergence criterion to use. The default value is 1e-5.
#'   \item \strong{power}: (Numeric) Raise factor loadings the the n-th power in the promax rotation. The default value is 4.
#'   \item \strong{norm}: (Logical) A logical value indicating whether to compute Kaiser normalization. Default is FALSE
#'   \item \strong{maxITR}: (Numeric) The maximum number of allowed iterations for convergence before stopping the rotation. The default value is 15,000 iterations.
#' }
#' @param faControl (List) A list of optional parameters passed to the factor extraction (\code{\link{faX}}) function.
#' \itemize{
#'   \item \strong{treatHeywood}: (Logical) In fals, if treatHeywood is true, a penalized least squares function is used to bound the communality estimates below 1.0. The default is TRUE.
#'   \item \strong{nStart}: (Numeric) In faml, determine the number of starting values to try. The default is 10 start values.
#'   \item \strong{maxCommunality}: (Numeric) In faml, set the maximum communality value for the estimated solution. The default maximum is .995.
#'   \item \strong{epsilon}: (Numeric) In fapa, the numeric threshold designating when the algorithm has converged. The default value is 1e-4.
#'   \item \strong{communality}: (Character) In fapa, the routine requires an initial communality estimate. Select how communalities are initially estimated. The default is squared multiple correlation ("SMC").
#'   \itemize{
#'     \item \strong{"SMC"}: Initial communalities are estimated by taking the squared multiple correlations of each indicator after regressing the indicator on the remaining variables. The following equation is employed to find the squared multiple correlation: \eqn{1 - 1 / diag(R^-1)}.
#'     \item \strong{"maxRsqr"}: Initial communalities equal the largest squared correlation in each column of the correlation matrix.
#'     \item \strong{"unity"}: Initial communalities equal 1.0 for all variables.
#'   }
#'   \item \strong{maxITR}: (Numeric) In fapa, the maximum number of iterations to reach convergence. The default is 15,000
#' }
#'
#' @return This function iterates the Schmid-Leiman target rotation and returns several relevant output.
#' \itemize{
#'   \item \strong{loadings}: (Matrix) The bifactor solution obtain from the SLi procedure.
#'   \item \strong{iterations}: (Numeric) The number of iterations required for convergence
#'   \item \strong{rotateControl}: (List) A list of the control parameters passed to the \code{\link{rotate}} function.
#'   \item \strong{faControl}: (List) A list of optional parameters passed to the factor extraction (\code{\link{faX}}) function.
#' }
#'
#' @author
#' \itemize{
#'   \item Casey Giordano (Giord023@umn.edu)
#'   \item Niels G. Waller (nwaller@umn.edu)
#'}
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
           urLambda      = NULL,      ## Unrotated factor structure
           freelyEstG    = TRUE,      ## Freely estimate all general factor loadings?
           gFac          = 1,         ## Column(s) of general factor(s)
           numberStarts  = 10,        ## num times to find minimum
           maxITR        = 20,        ## Max num of iterations
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


    ## Must specify either R or urLambda
    if ( is.null(urLambda) & is.null(R) ) {
      stop("Must specify either a correlation matrix (for factor extraction) or an unrotated factor structure.")
    } # END if ( is.null(urLambda) & is.null(R) )

    ## Check rotateControl

    ## rotateControl specified?
    if ( is.null(rotateControl) ) {

      ## Set the default rotateControl values (before user-specification of values)
      cnDefault <- list(gamma   = 0,
                        delta   = .01,
                        kappa   = 0,
                        k       = NULL,
                        epsilon = 1e-5,
                        power   = 4,
                        norm    = FALSE,
                        maxITR  = 15000)

      ## Used as func output, if it is NULL, need to specify
      ## Full error checking takes place in rotate()
      rotateControl <- cnDefault
    } # END if ( is.null(rotateControl) )

    ## Check faControl

    ## If faControl is not specified, give it the defaults (used for func output)
    if ( is.null(faControl) ) {

      ## Set the default values of all control agruments
      cnFA <- list(treatHeywood   = TRUE,
                   nStart         = 10,
                   maxCommunality = .995,
                   epsilon        = 1e-4,
                   communality    = "SMC",
                   maxITR         = 15000)

      ## Used as func output, if it is NULL, need to specify
      ## Full error checking takes place in faX()
      faControl <- cnFA

    } # END if ( is.null(faControl) )

    ## If digits is not specified, give it an arbitrarily large value
    if ( is.null(digits) ) {
      digits <- 100
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
    if ( is.null(urLambda) ) {

      ## Extract a full-rank factor structure matrix
      faOut <-
        faX(R          = R,
            numFactors = numberFactors,
            facMethod  = facMethod,
            faControl  = faControl)

      ## Save the unrotated factor loadings
      fact.loadings <- faOut$loadings[]

      ## Check for Heywood case
      if ( faOut$Heywood == TRUE ) {
        warning("There is a Heywood case in the unrotated, full-rank factor structure matrix.")
      } # END if ( faOut$Heywood == TRUE )

    } else {

      ## Use the pre-specified unrotated factor structure
      fact.loadings <- urLambda

    } # END if (is.null(urLambda)) {

    ## Begin a count for the number of iterations used for convergence
    numIters <- 0

    ## Iterate until convergence OR hitting maxITR limit
    for (iter in 1:maxITR){

      ## Procrustes rotation using the targetOld patterned matrix
      SLt <-   rotate(lambda        = fact.loadings,
                      rotate      = "targetT",
                      targetMatrix  = targetOld, ## This is updated each loop
                      numberStarts  = numberStarts,
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
    if (iter == maxITR & iter > 1) {

      warning("The solution did not converge in the maximum allotted iterations.
              \n Try increasing the maxITR argument.")

    } # END if (iter == maxITR & iter > 1) {

    ## Extract the final estimated solution
    SLi.loadings <- SLt$loadings

    ## Final output list
    ##    [[1]]: the final iterated Schmid-Leiman matrix of loadings
    ##    [[2]]: number of iterations to achieve final output

    list(loadings      = round(SLi.loadings, digits),
         iterations    = round(numIters, digits),
         rotateControl = rotateControl,
         faControl     = faControl)

  } # END SchmidLeimanIterated()
