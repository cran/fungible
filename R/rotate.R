#' Automatic Factor Rotation from Random Starting Configurations
#'
#' This function conducts factor rotations from a user-specified number of random (orthogonal) starting configurations. Based on the resulting discrepancy function, the function will determine the number of local minima and, among these local solutions, will find the "global minimum" (i.e., the minimal discrepancy value from a finite number of solutions).
#'
#' @param R (Matrix) A correlation matrix.
#' @param numFactors (Numeric) The number of factors to extract for subsequent rotation.
#' @param facMethod (Character) The method used for factor extraction. The supported options are "fals" for unweighted least squares, "faml" for maximum likelihood, and "fapa" for iterated principal axis factoring. The default method is "fals".
#' \itemize{
#'   \item \strong{"fals"}: Factors are extracted using the unweighted least squares estimation procedure using the \code{\link[fungible]{fals}} function.
#'   \item \strong{"faml"}: Factors are extracted using the maximum likelihood estimation procedure using the \code{\link[stats]{factanal}} function.
#'   \item \strong{"fapa"}: Factors are extracted using the iterated principal axis factoring estimation procedure using the \code{\link{fapa}} function.
#' }
#' @param lambda (Matrix) An unrotated factor-structure matrix to be rotated.
#' @param rotate (Character) Designate which rotation algorithm to apply. The following are available rotation options: "oblimin", "quartimin", "targetT", "targetQ", "oblimax", "entropy", "quartimax", "varimax", "simplimax", "bentlerT", "bentlerQ", "tandemI", "tandemII", "geominT", "geominQ", "cfT", "cfQ", "infomaxT", "infomaxQ", "mccammon", "bifactorT", "bifactorQ".
#' @param targetMatrix (Matrix) If a target rotation is desired, provide the target to be rotated toward. Note that a matrix of numeric values will yield a fully-specified Procrustes solution. To conduct Browne's partially-specified Procrustes rotation, ensure that all factor loadings to be freely estimated are designated by "NA" values. All non-salient values are typically specified as 0 (zero).
#' @param numberStarts (Numeric) The number of random (orthogonal) starting configurations to use. The default value is 10 random starts
#' @param digits (Numeric) The number of digits to round all output to. The default is no rounding.
#' @param keepAll (Logical) A logical value indicating whether to keep the solutions for each random starting configuration (i.e., all factor loadings, factor correlations, and their evaluated discrepancy function). The default value is TRUE.
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
#' @details
#' \itemize{
#'   \item \strong{Global Minimum}: This function uses several random starting configurations for factor rotations in an attempt to find the global minimum solution. However, this function is not guaranteed to find the global minimum. Furthermore, the global minimum solution need not be more psychologically interpretable than any of the local solutions (cf. Rozeboom, 1992). As is recommended, our function returns the purality of local solutions so users can make their own judgements.
#' }
#'
#' @return The rotate function will produce a lot of output in addition to the rotated factor pattern matrix and the factor correlations.
#' \itemize{
#' \item \strong{loadings}: (Matrix) The rotated factor solution with the lowest* evaluated discrepancy function. *This solution has the lowest discrepancy function \emph{of the examined random starting configurations}. It is not guarenteed to find the "true" global minimum. Note that multiple (or even all) local solutions can have the same discrepancy functions.
#' \item \strong{Phi}: (Matrix) The factor correlations of the rotated factor solution with the lowest* evaluated discrepancy function (see * in the loadings description above)
#' \item \strong{rotateControl}: (List) A list of the control parameters passed to the rotation (\code{\link{rotate}}) function.
#' \item \strong{faControl}: (List) A list of the control parameters passed to the factor extraction (\code{\link{faX}}) function.
#' \item \strong{localLoadings}: (List) A list of all of the local factor pattern matrices. If keepAll is FALSE, this will return NULL.
#' \item \strong{localPhi}: (List) A list of all of the local factor correlation matrices. If keepAll is FALSE, this will return NULL.
#' \item \strong{localDiscrepancy}: (List) A list of all of the minimized discrepancy functions. If keepAll is FALSE, this will return NULL.
#' \item \strong{minDiscrepancy}: (Scalar) A value, ranging from 1:numberStarts, indicating which of the local solutions has the smallest discrepancy function value. If keepAll is FALSE, this will return NULL.
#' }
#'
#'
#' @references Browne, M. W. (2001). An overview of analytic rotation in exploratory factor analysis. \emph{Multivariate Behavioral Research, 36}(1), 111-150.
#' @references Mansolf, M., & Reise, S. P. (2016). Exploratory bifactor analysis: The Schmid-Leiman orthogonalization and Jennrich-Bentler analytic rotations. \emph{Multivariate Behavioral Research, 51}(5), 698-717.
#' @references Rozeboom, W. W. (1992). The glory of suboptimal factor rotation: Why local minima in analytic optimization of simple structure are more blessing than curse. \emph{Multivariate Behavioral Research, 27}(4), 585-599.
#'
#' @author
#' \itemize{
#'   \item Casey Giordano (Giord023@umn.edu)
#'   \item Niels G. Waller (nwaller@umn.edu)
#'}
#'
#' @examples
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
#' ## Model-implied correlation (covariance) matrix
#' R <- lambda %*% Phi %*% t(lambda)
#' diag(R) <- 1
#'
#' ## 100 promax rotations from least squares factor extraction
#' Out1 <- rotate(R             = R,
#'                numFactors    = 3,
#'                facMethod     = "fals",
#'                rotate        = "promaxQ",
#'                numberStarts  = 100,
#'                keepAll       = TRUE,
#'                rotateControl = list(power = 4,
#'                                     norm  = TRUE),
#'                faControl     = list(treatHeywood = TRUE),
#'                digits        = 2)$loadings
#'
#' ## 10 oblimin rotations from maximum likelihood factor extraction
#' Out2 <- rotate(R             = R,
#'                numFactors    = 3,
#'                facMethod     = "faml",
#'                rotate        = "oblimin",
#'                numberStarts  = 10,
#'                keepAll       = FALSE,
#'                rotateControl = list(gamma   = 0,
#'                                     epsilon = 1e-6),
#'                faControl     = list(nStart         = 10,
#'                                     maxCommunality = .99),
#'                digits        = 2)$loadings
#'
#' @export


rotate <-
  function(R             = NULL,   ## Correlation matrix for extraction
           numFactors    = NULL,   ## Number of factors to extract
           facMethod     = "fals", ## Factor extraction method
           lambda        = NULL,   ## Unrotated loading matrix
           rotate        = NULL,   ## Which rotation to use
           targetMatrix,           ## If targetT, matrix of specified loadings
           numberStarts  = 10,     ## Number of attempts to try
           digits        = NULL,   ## Round all output
           keepAll       = TRUE,   ## Keep output from local solutions
           rotateControl = NULL,   ## Control rotation tuning parameters
           faControl     = NULL) { ## Control factor extraction parameters

    ## ~~~~~~~~~~~~~~~~~~ ##
    #### Error checking ####
    ## ~~~~~~~~~~~~~~~~~~ ##

    ## Check rotation

    ## That a rotation is specified
    if ( is.null(rotate) ) {
      stop("The 'rotate' argument must be specified. See ?rotate for available options")
    } # END if ( is.null(rotate) )

    ## Character string of all available rotate options
    PlausibleRotations <- c("oblimin",   "quartimin", "targetT",
                            "targetQ",   "oblimax",   "entropy",
                            "quartimax", "varimax",   "simplimax",
                            "bentlerT",  "bentlerQ",  "tandemI",
                            "tandemII",  "geominT",   "geominQ",
                            "promaxQ",   "cfT",       "cfQ",
                            "infomaxT",  "infomaxQ",  "mccammon",
                            "bifactorT", "bifactorQ")

    ## If rotate argument is incorrectly specified, return error
    if ( rotate %in% PlausibleRotations == FALSE ) {
      stop("The 'rotate' argument is incorrectly specified. Check for spelling errors (case sensitive).")
    } # END if (rotate %in% PlausibleRotations == FALSE) {

    ## Check R

    ## Make sure either corr mat or factor structure is supplied
    if ( is.null(R) & is.null(lambda) ) {
      stop("Must specify either 'R' or 'lambda' arguments.")
    } # END if ( is.null(R) & is.null(lambda) )

    ## If factor analysis is run, must specify numFactors
    if ( !is.null(R) & is.null(numFactors) ) {
      stop("A correlation matrix is supplied but the user needs to specify the number of factors to extract.")
    } # END if ( !is.null(R) & is.null(numFactors) )

    ## If a digits argument is not supplied, give an arbitrarily large one
    if ( is.null(digits) ) {
      digits <- 100
    } # END if ( is.null(digits) )

    ## Assign the default values for the rotateControl list
    cnDefault <- list(gamma   = 0,
                      delta   = .01,
                      kappa   = 0,
                      k       = nrow(lambda),
                      epsilon = 1e-5,
                      power   = 4,
                      norm    = FALSE,
                      maxITR  = 15000)

    ## Test that rotateControl arguments are correctly specified
    if ( !is.null(rotateControl) ) {

      ## Number of names in the default rotateControl list
      cnLength <- length( names(cnDefault) )

      ## Number of all unique names across user rotateControl and default rotateControl lists
      allLength <- length( unique( c( names(rotateControl), names(cnDefault) ) ) )

      ## If the lengths differ, it is because of spelling issue
      if (cnLength != allLength) {

        ## Find the incorrect arg/s
        incorrectArgs <- which( (names(rotateControl) %in% names(cnDefault) ) == FALSE)

        # stop("The following arguments are not valid inputs for the list of rotateControl arguments: ", paste0( paste(c(names(rotateControl)[incorrectArgs]), collapse = ", "), "." ) )
        stop(paste("The following arguments are not valid inputs for the list of rotateControl arguments:", paste(c(names(rotateControl)[incorrectArgs]), collapse = ", "), collapse = ": " ) )

      } # END if (cnLength != allLength)

    } # END if ( is.null(rotateControl) )


    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    #### Begin body of the function ####
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

    ## Change the default values based on user-specified rotateControl arguments
    cnDefault[names(rotateControl)] <- rotateControl

    ## Generate a random orthonormal starting matrix
    randStart <- function(dimension) {
      qr.Q(qr(matrix(rnorm(dimension^2), dimension, dimension )))
    } # END randStart <- function(dimension) {

    ## If lambda is not supplied, extract it
    if ( is.null(lambda) ) {

      ## Call faX function to extract factor structure matrix
      lambda <- faX(R          = R,
                    numFactors = numFactors,
                    facMethod  = facMethod,
                    faControl  = faControl)$loadings

    } # END if ( is.null(lambda) )

    ## Determine size of the random starting matrix
    matrixDim <- ncol(lambda)

    ## Pre-allocate a list for the different attempts
    starts <- vector("list",
                     numberStarts)

    ## For each list element, do the specified rotate and save all output
    starts <- lapply(starts, function(x)
      switch(rotate,
             "oblimin" = {
               GPArotation::oblimin(lambda,
                                    Tmat      = randStart(matrixDim),
                                    gam       = cnDefault$gamma,
                                    normalize = cnDefault$norm,
                                    eps       = cnDefault$epsilon,
                                    maxit     = cnDefault$maxITR)
             },
             "quartimin" = {
               GPArotation::quartimin(lambda,
                                      Tmat      = randStart(matrixDim),
                                      maxit     = cnDefault$maxITR,
                                      eps       = cnDefault$epsilon,
                                      normalize = cnDefault$norm)
             },
             "targetT" = {
               GPArotation::targetT(lambda,
                                    Tmat      = randStart(matrixDim),
                                    Target    = targetMatrix,
                                    normalize = cnDefault$norm,
                                    eps       = cnDefault$epsilon,
                                    maxit     = cnDefault$maxITR)
             },
             "targetQ" = {
               GPArotation::targetQ(lambda,
                                    Tmat      = randStart(matrixDim),
                                    Target    = targetMatrix,
                                    normalize = cnDefault$norm,
                                    eps       = cnDefault$epsilon,
                                    maxit     = cnDefault$maxITR)
             },
             "oblimax" = {
               GPArotation::oblimax(lambda,
                                    Tmat      = randStart(matrixDim),
                                    normalize = cnDefault$norm,
                                    eps       = cnDefault$epsilon,
                                    maxit     = cnDefault$maxITR)
             },
             "entropy" = {
               GPArotation::entropy(lambda,
                                    Tmat      = randStart(matrixDim),
                                    normalize = cnDefault$norm,
                                    eps       = cnDefault$epsilon,
                                    maxit     = cnDefault$maxITR)
             },
             "quartimax" = {
               GPArotation::quartimax(lambda,
                                      Tmat      = randStart(matrixDim),
                                      normalize = cnDefault$norm,
                                      eps       = cnDefault$epsilon,
                                      maxit     = cnDefault$maxITR)
             },
             "varimax" = {
               GPArotation::Varimax(lambda,
                                    Tmat      = randStart(matrixDim),
                                    normalize = cnDefault$norm,
                                    eps       = cnDefault$epsilon,
                                    maxit     = cnDefault$maxITR)
             },
             "simplimax" = {
               GPArotation::simplimax(lambda,
                                      Tmat      = randStart(matrixDim),
                                      normalize = cnDefault$norm,
                                      eps       = cnDefault$epsilon,
                                      k         = cnDefault$k,
                                      maxit     = cnDefault$maxITR)
             },
             "bentlerT" = {
               GPArotation::bentlerT(lambda,
                                     Tmat      = randStart(matrixDim),
                                     maxit     = cnDefault$maxITR,
                                     eps       = cnDefault$epsilon,
                                     normalize = cnDefault$norm)
             },
             "bentlerQ" = {
               GPArotation::bentlerQ(lambda,
                                     Tmat      = randStart(matrixDim),
                                     maxit     = cnDefault$maxITR,
                                     eps       = cnDefault$epsilon,
                                     normalize = cnDefault$norm)
             },
             "tandemI" = {
               GPArotation::tandemI(lambda,
                                    Tmat      = randStart(matrixDim),
                                    maxit     = cnDefault$maxITR,
                                    eps       = cnDefault$epsilon,
                                    normalize = cnDefault$norm)
             },
             "tandemII" = {
               GPArotation::tandemII(lambda,
                                     Tmat      = randStart(matrixDim),
                                     maxit     = cnDefault$maxITR,
                                     eps       = cnDefault$epsilon,
                                     normalize = cnDefault$norm)
             },
             "geominT" = {
               GPArotation::geominT(lambda,
                                    Tmat      = randStart(matrixDim),
                                    delta     = cnDefault$delta,
                                    normalize = cnDefault$norm,
                                    eps       = cnDefault$epsilon,
                                    maxit     = cnDefault$maxITR)
             },
             "geominQ" = {
               GPArotation::geominQ(lambda,
                                    Tmat      = randStart(matrixDim),
                                    delta     = cnDefault$delta,
                                    normalize = cnDefault$norm,
                                    eps       = cnDefault$epsilon,
                                    maxit     = cnDefault$maxITR)
             },
             "promaxQ" = {
               promaxQ(lambda  = lambda,
                       power   = cnDefault$power,
                       norm    = cnDefault$norm,
                       epsilon = cnDefault$epsilon,
                       maxITR  = cnDefault$maxITR)
             },
             "cfT" = {
               GPArotation::cfT(lambda,
                                Tmat      = randStart(matrixDim),
                                kappa     = cnDefault$kappa,
                                maxit     = cnDefault$maxITR,
                                eps       = cnDefault$epsilon,
                                normalize = cnDefault$norm)
             },
             "cfQ" = {
               GPArotation::cfQ(lambda,
                                Tmat      = randStart(matrixDim),
                                kappa     = cnDefault$kappa,
                                eps       = cnDefault$epsilon,
                                normalize = cnDefault$norm,
                                maxit     = cnDefault$maxITR)
             },
             "infomaxT" = {
               GPArotation::infomaxT(lambda,
                                     Tmat      = randStart(matrixDim),
                                     normalize = cnDefault$norm,
                                     eps       = cnDefault$epsilon,
                                     maxit     = cnDefault$maxITR)
             },
             "infomaxQ" = {
               GPArotation::infomaxQ(lambda,
                                     Tmat      = randStart(matrixDim),
                                     normalize = cnDefault$norm,
                                     eps       = cnDefault$epsilon,
                                     maxit     = cnDefault$maxITR)
             },
             "mccammon" = {
               GPArotation::mccammon(lambda,
                                     Tmat      = randStart(matrixDim),
                                     normalize = cnDefault$norm,
                                     eps       = cnDefault$epsilon,
                                     maxit     = cnDefault$maxITR)
             },
             "bifactorT" = {
               GPArotation::bifactorT(lambda,
                                      Tmat      = randStart(matrixDim),
                                      normalize = cnDefault$norm,
                                      eps       = cnDefault$epsilon,
                                      maxit     = cnDefault$maxITR)
             },
             "bifactorQ" = {
               GPArotation::bifactorQ(lambda,
                                      Tmat      = randStart(matrixDim),
                                      normalize = cnDefault$norm,
                                      eps       = cnDefault$epsilon,
                                      maxit     = cnDefault$maxITR)
             }))

    ## Find the minimum value of each rotate, save it
    min.val <- rep(NA, numberStarts)

    ## Save the minimized discrepancy function for each attempt
    ##    $Table[,2] is the value of criterion at each iteration, which.min
    ##    finds the minimum criterion value across each attempt

    ## Evaluate the minimum disc functions to find smallest value
    if (rotate == "promaxQ") {

      ## For each random start, find the evaluated discrepancy function
      DiscrepFunc <- lapply(starts, function(attempt) attempt$vmaxDiscrepancy)

    } else {

      ## For each random start, find the evaluated discrepancy function
      DiscrepFunc <- lapply(starts, function(attempt) min(attempt$Table[, 2]))

    } # END if (rotate == "promaxQ")

    ## Of all random configs, determinine which has the lowest criterion value
    minDiscrepFunc <- which.min(DiscrepFunc)

    ## Select the loading matrix with the global minimum
    rotated.lambda <- starts[[minDiscrepFunc]]$loadings
    rotated.phi    <- starts[[minDiscrepFunc]]$Phi

    ## List output:
    ##    [[loadings]]: the rotate factor solution. Of the random starting configurations, it is the rotation with the lowest evaluated discrepancy function. It is possible that multiple (or even all) local solutions are exactly equal.
    ##    [[Phi]]: rotated phi matrix that had minimum discrepancy function
    ##    [[3]]: vector of logicals stating whether each of numStarts attempts converged
    ##    [[4]]: logical stating whether minimum of all attempted rotations converged

    rotation.output <-
      list(loadings        = round(rotated.lambda, digits),
           Phi             = rotated.phi,
           rotateControl   = cnDefault,
           faControl       = faControl,
           localLoadings   = lapply(starts, function(x) x$loadings),
           localPhi        = lapply(starts, function(x) x$Phi),
           localDiscepancy = DiscrepFunc,
           minDiscrepancy  = minDiscrepFunc)

    ## Round Phi manually (sometimes it is NULL, don't round it then)
    if ( is.matrix(rotation.output$Phi) ) {
      round(rotation.output$Phi, 2)
    } # END if ( is.matrix(rotation.output$Phi) )

    ## Turn values to NULL to reduce the length of the output
    if ( keepAll == FALSE ) {

      ## Turn the following lists to a NULL value
      rotation.output$localLoadings <-
        rotation.output$localPhi <-
        rotation.output$localDiscepancy <- NULL
    } # END if ( keepAll == FALSE )

    ## Return the final output
    rotation.output

  } ## END rotate()
