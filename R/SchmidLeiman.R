#' Schmid-Leiman Orthogonalization to a (Rank-Deficient) Bifactor Structure
#'
#' The Schmid-Leiman (SL) procedure orthogonalizes a higher-order factor structure into a rank-deficient bifactor structure. The Schmid-Leiman method is a generalization of Thomson's orthogonalization routine.
#'
#' @param R (Matrix) A correlation matrix.
#' @param numFactors (Vector) The number of latent factors at each level of analysis. For example, c(3, 1) estimates three latent factors in the first-order common factor model and one latent factor in the second-order common factor model (i.e., 3 group factors and 1 general factor). This function can orthogonalize up to (and including) a three-order factor solution.
#' @param facMethod (Character) The method used for factor extraction (see \code{\link{faX}} for more details). The supported options are "fals" for unweighted least squares, "faml" for maximum likelihood, and "fapa" for iterated principal axis factoring. The default method is "fals".
#' \itemize{
#'   \item \strong{"fals"}: Factors are extracted using the unweighted least squares estimation procedure using the \code{\link[fungible]{fals}} function.
#'   \item \strong{"faml"}: Factors are extracted using the maximum likelihood estimation procedure using the \code{\link[stats]{factanal}} function.
#'   \item \strong{"fapa"}: Factors are extracted using the iterated principal axis factoring estimation procedure using the \code{\link{fapa}} function.
#' }
#' @param rotate (Character) Designate which rotation algorithm to apply. See the \code{\link{rotate}} function for more details about possible rotations. A geomin rotation is the default.
#' @param numberStarts (Numeric) Designate the number of random starting configurations for each factor rotation.
#' @param rescaleH2 (Numeric) If a Heywood case is detected at any level of the higher-order factor analyses, rescale the communality value to continue with the matrix algebra. When a Heywood case occurs, the uniquenesses (i.e., specific-factor variances) will be negative and the SL orthogonalization of the group factors is no longer correct.
#' @param digits (Numberic) Number of digits to round factor loadings in the resulting bifactor solution.
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
#' @return
#' \itemize{
#'   \item \strong{L1}: (Matrix) The first-order (oblique) factor pattern matrix.
#'   \item \strong{L2}: (Matrix) The second-order (oblique) factor pattern matrix.
#'   \item \strong{L3}: (Matrix, NULL) The third-order (oblique) factor pattern matrix (if applicable).
#'   \item \strong{Phi1}: (Matrix) The first-order factor correlation matrix.
#'   \item \strong{Phi2}: (Matrix) The second-order factor correlation matrix.
#'   \item \strong{Phi3}: (Matrix, NULL) The third-order factor pattern matrix (if applicable).
#'   \item \strong{Usq1}: (Matrix) The first-order factor uniquenesses (variances).
#'   \item \strong{Usq2}: (Matrix) The second-order factor uniquenesses (variances).
#'   \item \strong{Usq3}: (Matrix, NULL) The third-order factor uniquenesses (variances) (if applicable).
#'   \item \strong{B}: (Matrix) The resulting Schmid-Leiman transformation.
#'   \item \strong{rotateControl}: (List) A list of the control parameters passed to the \code{\link{rotate}} function.
#'   \item \strong{faControl}: (List) A list of optional parameters passed to the factor extraction (\code{\link{faX}}) function.
#'}
#'
#' @references Abad, F. J., Garcia-Garzon, E., Garrido, L. E., & Barrada, J. R. (2017). Iteration of partially specified target matrices: application to the bi-factor case. \emph{Multivariate Behavioral Research, 52}(4), 416-429.
#' @references Giordano, C. & Waller, N. G. (under review). Recovering bifactor models: A comparison of seven methods.
#' @references Schmid, J., & Leiman, J. M. (1957). The development of hierarchical factor solutions. \emph{Psychometrika, 22}(1), 53-61.
#'
#' @author
#' \itemize{
#'   \item Casey Giordano (Giord023@umn.edu)
#'   \item Niels G. Waller (nwaller@umn.edu)
#'}
#'
#' @examples
#' ## Dataset used in Schmid & Leiman (1957) rounded to 2 decimal places
#' SLdata <-
#'   matrix(c(1.0, .72, .31, .27, .10, .05, .13, .04, .29, .16, .06, .08,
#'            .72, 1.0, .35, .30, .11, .06, .15, .04, .33, .18, .07, .08,
#'            .31, .35, 1.0, .42, .08, .04, .10, .03, .22, .12, .05, .06,
#'            .27, .30, .42, 1.0, .06, .03, .08, .02, .19, .11, .04, .05,
#'            .10, .11, .08, .06, 1.0, .32, .13, .04, .11, .06, .02, .03,
#'            .05, .06, .04, .03, .32, 1.0, .07, .02, .05, .03, .01, .01,
#'            .13, .15, .10, .08, .13, .07, 1.0, .14, .14, .08, .03, .04,
#'            .04, .04, .03, .02, .04, .02, .14, 1.0, .04, .02, .01, .01,
#'            .29, .33, .22, .19, .11, .05, .14, .04, 1.0, .45, .15, .17,
#'            .16, .18, .12, .11, .06, .03, .08, .02, .45, 1.0, .08, .09,
#'            .06, .07, .05, .04, .02, .01, .03, .01, .15, .08, 1.0, .42,
#'            .08, .08, .06, .05, .03, .01, .04, .01, .17, .09, .42, 1.0),
#'          nrow = 12, ncol = 12, byrow = TRUE)
#'
#' Out1 <- SchmidLeiman(R          = SLdata,
#'                      numFactors = c(6, 3, 1),
#'                      digits     = 2)$B
#'
#' ## An orthogonalization of a two-order structure
#' bifactor <- matrix(c(.16, .37, .00, .00,
#'                      .18, .41, .00, .00,
#'                      .21, .48, .00, .00,
#'                      .36, .00, .55, .00,
#'                      .21, .00, .32, .00,
#'                      .36, .00, .55, .00,
#'                      .47, .00, .00, .48,
#'                      .40, .00, .00, .40,
#'                      .39, .00, .00, .39),
#'                    nrow = 9, ncol = 3, byrow = TRUE)
#'
#' ## Model-implied correlation (covariance) matrix
#' R <- bifactor %*% t(bifactor)
#'
#' ## Unit diagonal elements
#' diag(R) <- 1
#'
#' Out1 <- SchmidLeiman(R          = R,
#'                      numFactors = c(3, 1),
#'                      digits     = 2)$B
#'
#' @export


SchmidLeiman <-
  function(R,                         ## Correlation matrix
           numFactors,                ## Number of grp factors per level
           facMethod     = "fals",    ## Factor extraction method
           rotate        = "geominQ", ## Oblique rotation
           numberStarts  = 10,        ## Number of random starts
           rescaleH2     = .98,       ## Rescale communalities
           digits        = NULL,      ## Round all output
           rotateControl = NULL,      ## Arguments passed to rotate
           faControl     = NULL){     ## Arguments passed to faX

    ## ~~~~~~~~~~~~~~~~~~~~~~ ##
    #### Internal Functions ####
    ## ~~~~~~~~~~~~~~~~~~~~~~ ##

    ## Rescale the factor loadings when a Heywood case occurs
    commRescale <- function(lambda,
                            threshold = 1.0,
                            rescaleTo = .98) {
      ## Rescale the factor structure matrix when Heywood case is detected
      hsq.new <-
        hsq.old <-
        apply(lambda^2, 1, sum)

      ## communalities that exceed threshold will be rescaled
      hsq.new[which(hsq.old >= threshold)] <- rescaleTo

      ## Norm lambda by dividing each row by its communality
      ## Then pre-multiply by the new (sqrt) communality diagonal matrix
      newLambda <- diag(sqrt(hsq.new)) %*%
        sqrt(diag(1 / hsq.old)) %*%
        lambda

      newLambda
    } # END commRescale


    ## ~~~~~~~~~~~~~~~~~~ ##
    #### Error Checking ####
    ## ~~~~~~~~~~~~~~~~~~ ##

    ## Is correlation matrix symmetric
    if ( !isSymmetric(R) ) {
      stop("The user-defined correlation matrix is not symmetric")
    } # END if (!isSymmetric(R)) {

    ## Number of factors correctly specified?
    if (length(numFactors) <= 1) {
      stop("The 'numFactors' argument must be a vector of 2 or 3 values.")
    } # END if (length(numFactors) <=1)

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

    ## If digits is not specified, give arbitrarily large value
    if ( is.null(digits) ) {
      digits <- 100
    } # END if ( is.null(digits) )

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ##         First-level hierarchical factor analysis         ##
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

    ## Extract 1st level of the hierarchy

    ## Extract the factor solution
    L1Out <- faX(R          = R,
                         numFactors = numFactors[1],
                         facMethod  = facMethod,
                         faControl  = faControl)

    ## Save the factor loadings for rotation
    L1 <- L1Out$loadings

    ## Rotate the factor structure
    firstRot <- rotate(lambda        = L1,
                         rotate      = rotate,
                         numberStarts  = numberStarts,
                         rotateControl = rotateControl)

    ## Store the Phi matrix and rotated loadings matrix
    L1 <- firstRot$loadings   ## Factor pattern
    R1 <- firstRot$Phi        ## Factor correlations

    ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## Check for Heywood cases ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~ ##

    ## Save level-one communality values
    HSquare <- L1Out$h2

    ## If any communalities are greater than lambda > 1.0, re-scale problem items
    if ( L1Out$Heywood == TRUE ) {

      ## Rescale the communalities larger than 'threshold' to 'rescaleTo'
      L1 <- commRescale(lambda    = L1,
                        threshold = 1.0,
                        rescaleTo = rescaleH2)

      message("A Heywood case was detected in the first-order factor pattern. Communalities have been rescaled.")

      ## Find new communality values after rescaling
      HSquare <- apply(L1^2, 1, sum)

    } # END if ( L1Out$Heywood == TRUE )

    ## Lower bound of uniqueness is defined at .05 (1 - newCommunality)
    U2_1 <- diag(1 - as.vector(HSquare))

    ## Save sqrt of U2_1 for later orthogonalization
    U_1 <- sqrt(U2_1)

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ##         Second-level hierarchical factor analysis        ##
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

    ## Unrotated factor loadings of the 2nd level of analysis
    L2Out <- faX(R          = R1,
                         numFactors = numFactors[2],
                         facMethod  = facMethod,
                         faControl  = faControl)

    ## Save the factor loadings
    L2 <- L2Out$loadings[]

    ## Only needs to rotate obliquely if more than 1 factor is extracted
    if (numFactors[2] > 1) {

      ## Rotate the factor structure
      secondRot <- rotate(lambda        = L2,
                            rotate      = rotate,
                            numberStarts  = numberStarts,
                            rotateControl = rotateControl)

      ## Store the Phi matrix and rotated loadings matrix
      L2 <- secondRot$loadings

      ## Detect if Heywood cases are present
      HSquare.2 <- L2Out$h2

      ## Stop if a Heywood case is detected
      if ( L2Out$Heywood == TRUE ) {

        ## Rescale the communalities larger than 'threshold' to 'rescaleTo'
        L2 <- commRescale(lambda    = L2,
                          threshold = 1.0,
                          rescaleTo = rescaleH2)

        message("A Heywood case was detected in the second-order factor pattern. Communalities have been rescaled.")

        ## Find rescaled communalities
        HSquare.2 <- apply(L2^2, 1, sum)

      } # END if ( L2Out$Heywood == TRUE )

      ## Factor uniqueness values
      U2_2 <- diag(1 - as.vector(HSquare.2))

      ## Save factor standard deviations for orthogonalization
      U_2 <- sqrt(U2_2)

      ## Save the phi matrix
      R2 <- secondRot$Phi

    } else {  #END if(numFactors[2] > 1)

      ## Check for Heywood cases
      HSquare.2 <- L2Out$h2

      if ( L2Out$Heywood == TRUE ) {
        ## Rescale the communalities larger than 'threshold' to 'rescaleTo'
        L2 <- commRescale(lambda    = L2,
                          threshold = 1.0,
                          rescaleTo = rescaleH2)

        message("A Heywood case was detected in the second-order factor pattern.
              Communalities have been rescaled.")

        ## Find rescaled communality
        HSquare.2 <- L2^2

      } # END if ( any(HSquare >= .99) )

      ## Factor uniqueness values
      U2_2 <- diag(1 - as.vector(HSquare.2))

      ## No phi matrix when only 1 factor is present, set at unity
      R2 <- 1

      ## Save for orthogonalization later
      U_2 <- sqrt(U2_2)

    } # END if (numFactors[2] > 1)

    ## Orthogonalize, Final step in Schmid & Leiman ('57) with 2 levels
    B.out <- cbind(L1 %*% L2, L1 %*% U_2)



    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ##         Third-level hierarchical factor analysis         ##
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

    ## If a 3rd level of hierarchy exists, begin the process

    if (length(numFactors) == 3) {

      ## Unrotated factor analysis of the 3rd order factor model
      L3Out <- faX(R          = R2,
                           numFactors = numFactors[3],
                           facMethod  = facMethod,
                           faControl  = faControl)

      ## Save the factor loadings
      L3 <- L3Out$loadings[]

      ## If multiple factors exist in the 3rd level, rotate obliquely
      if (numFactors[3] > 1) {

        ## Conduct oblique rotation if multiple 3rd-order factors exist
        thirdRot <- rotate(lambda        = L3,
                             rotate      = rotate,
                             numberStarts  = numberStarts,
                             rotateControl = rotateControl)

        ## Obliquely rotated factor loadings in the 3rd level
        L3 <- thirdRot$loadings

        ## Compute the communality estimate
        HSquare.3 <- apply(L3^2, 1, sum)

        ## Detect if Heywood cases are present
        if ( any(HSquare.3 > 1.0) ) {

          ## Rescale the communalities larger than 'threshold' to 'rescaleTo'
          L3 <- commRescale(lambda    = L3,
                            threshold = 1.0,
                            rescaleTo = rescaleH2)

          message("A Heywood case was detected in the third-order factor pattern.
                Communalities have been rescaled.")

          ## Find rescaled communality
          HSquare.3 <- apply(L3^2, 1, sum)

        } # END if ( any(HSquare.3 > 1.0) )

        ## Faster to compute from saved output
        U2_3 <- diag(1 - as.vector(HSquare.3))

        ## Set the 3rd level phi matrix
        R3 <- thirdRot$Phi

      } else { ## if(numFactors[3] == 1)

        ## Compute the uniqueness matrix
        U2_3 <- diag(1 - as.vector(L3^2))

        ## Set phi matrix to unity (only 1 factor)
        R3 <- 1
      } # END if(numFactors[3] == 1)

      ## No matter how U2_3 is found, save the sqrt of it
      U_3 <- sqrt(U2_3)

      ## Start orthogonalization for 3 levels
      B3    <- cbind(L3, U_3)
      B2    <- cbind((L2 %*% B3), U_2)
      B.out <- L1 %*% B2

      ## Give the 3-level analysis some labels
      colnames(B.out) <-
        c(paste(rep("g", numFactors[3]), c(1:numFactors[3]), sep = ""),
          paste(rep("F", numFactors[2]), c(1:numFactors[2]), sep = ""),
          paste(rep("f", numFactors[1]), c(1:numFactors[1]), sep = ""))

      ## Names of items resemble that of initial correlation matrix <-
      rownames(B.out) <- rownames(R)
      # END if(length(numFactors) == 3))

    } else { ## Do this if only 2 levels exist

      ## Give the 2-level analysis some labels, ARBITRARY order
      colnames(B.out) <-
        c(paste(rep("g", numFactors[2]), c(1:numFactors[2]), sep = ""),
          paste(rep("F", numFactors[1]), c(1:numFactors[1]), sep = ""))

      ## Names of items resemble that of initial correlation matrix
      rownames(B.out) <- rownames(R)

    }## END if(length(numFactors) == 3){

    ## Available output
    SLOutput <-
      list(L1            = round(L1, digits),    ## level 1 loadings
           L2            = round(L2, digits),    ## level 2 loadings
           L3            = NULL,  ## level 3 loadings
           Phi1          = round(R1, digits),    ## level 1 phi matrix
           Phi2          = round(R2, digits),    ## level 2 phi matrix
           Phi3          = NULL,  ## level 3 phi matrix
           Usq1          = round(U2_1, digits),  ## level 1 uniquenesses
           Usq2          = round(U2_2, digits),  ## level 2 uniquenesses
           Usq3          = NULL,  ## level 3 uniquenesses
           B             = round(B.out, digits),
           rotateControl = rotateControl,
           faControl     = faControl) ## final bifactor solution

    ## If a 3-rd order solution is produced, populate the output of the 3rd level
    if ( length(numFactors) == 3 ) {

      SLOutput$L3   <- round(L3,   digits)
      SLOutput$Phi3 <- round(R3,   digits)
      SLOutput$Usq3 <- round(U2_3, digits)

    } # END if ( length(numFactors) == 3 )

    ## Return the list of output
    SLOutput

  }## END SchmidLeiman()
