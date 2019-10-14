#' Schmid-Leiman Orthogonalization to a (Rank-Deficient) Bifactor Structure
#'
#' The Schmid-Leiman (SL) procedure orthogonalizes a higher-order factor 
#' structure into a rank-deficient bifactor structure. The Schmid-Leiman method 
#' is a generalization of Thomson's orthogonalization routine.
#'
#' @param R (Matrix) A correlation matrix.
#' @param numFactors (Vector) The number of latent factors at each level of 
#' analysis. For example, c(3, 1) estimates three latent factors in the first-order 
#' common factor model and one latent factor in the second-order common factor 
#' model (i.e., 3 group factors and 1 general factor). This function can 
#' orthogonalize up to (and including) a three-order factor solution.
#' @param rotate (Character) Designate which rotation algorithm to apply. See 
#' the \code{\link{faMain}} function for more details about possible rotations. 
#' Defaults to rotate = "oblimin".
#' @param rescaleH2 (Numeric) If a Heywood case is detected at any level of the 
#' higher-order factor analyses, rescale the communality value to continue with 
#' the matrix algebra. When a Heywood case occurs, the uniquenesses (i.e., 
#' specific-factor variances) will be negative and the SL orthogonalization of 
#' the group factors is no longer correct.
#' @inheritParams faMain
#'
#' @details The obtained Schmid-Leiman (SL) factor structure matrix is rescaled if 
#'  its communalities differ from those of the  original first-order 
#'  solution (due to the presence of one or more Heywood cases in a solution of any order). 
#'  Rescaling will produce SL communalities that match those of 
#'  the original first-order  solution.
#'
#' @return
#' \itemize{
#'   \item \strong{L1}: (Matrix) The first-order (oblique) factor pattern matrix.
#'   \item \strong{L2}: (Matrix) The second-order (oblique) factor pattern matrix.
#'   \item \strong{L3}: (Matrix, NULL) The third-order (oblique) factor pattern 
#'   matrix (if applicable).
#'   \item \strong{Phi1}: (Matrix) The first-order factor correlation matrix.
#'   \item \strong{Phi2}: (Matrix) The second-order factor correlation matrix.
#'   \item \strong{Phi3}: (Matrix, NULL) The third-order factor pattern matrix 
#'   (if applicable).
#'   \item \strong{U1}: (Matrix) The square root of the first-order factor 
#'   uniquenesses (i.e., factor standard deviations).
#'   \item \strong{U2}: (Matrix) The square root of the second-order factor 
#'   uniquenesses (i.e., factor standard deviations).
#'   \item \strong{U3}: (Matrix, NULL) The square root of the third-order factor 
#'   uniquenesses (i.e., factor standard deviations) (if applicable).
#'   \item \strong{B}: (Matrix) The resulting Schmid-Leiman transformation.
#'   \item \strong{rotateControl}: (List) A list of the control parameters 
#'   passed to the \code{\link{faMain}} function.
#'   \item \strong{faControl}: (List) A list of optional parameters passed to 
#'   the factor extraction (\code{\link{faX}}) function.
#'   \item{strong{HeywoodFlag}}(Integer) An integer indicating whether one or more 
#'    Heywood cases were encountered during estimation.
#'}
#'
#' @references Abad, F. J., Garcia-Garzon, E., Garrido, L. E., & Barrada, J. R. 
#' (2017). Iteration of partially specified target matrices: application to the 
#' bi-factor case. \emph{Multivariate Behavioral Research, 52}(4), 416-429.
#' @references Giordano, C. & Waller, N. G. (under review). Recovering bifactor 
#' models: A comparison of seven methods.
#' @references Schmid, J., & Leiman, J. M. (1957). The development of 
#' hierarchical factor solutions. \emph{Psychometrika, 22}(1), 53-61.
#'
#' @family Factor Analysis Routines
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
#' bifactor <- matrix(c(.46, .57, .00, .00,
#'                      .48, .61, .00, .00,
#'                      .61, .58, .00, .00,
#'                      .46, .00, .55, .00,
#'                      .51, .00, .62, .00,
#'                      .46, .00, .55, .00,
#'                      .47, .00, .00, .48,
#'                      .50, .00, .00, .50,
#'                      .49, .00, .00, .49),
#'                    nrow = 9, ncol = 3, byrow = TRUE)
#'
#' ## Model-implied correlation (covariance) matrix
#' R <- bifactor %*% t(bifactor)
#'
#' ## Unit diagonal elements
#' diag(R) <- 1
#'
#' Out2 <- SchmidLeiman(R          = R,
#'                      numFactors = c(3, 1),
#'                      rotate     = "oblimin",
#'                      digits     = 2)$B
#'
#' @export

SchmidLeiman <-
  function(R,                         ## Correlation matrix
           numFactors,                ## Number of grp factors per level
           facMethod     = "fals",    ## Factor extraction method
           rotate        = "oblimin", ## Oblique rotation
           rescaleH2     = .98,       ## Rescale communalities
           digits        = NULL,      ## Round all output
           faControl     = NULL,      ## Arguments passed to faX
           rotateControl = NULL){     ## Arguments passed to rotate
    
    ## ~~~~~~~~~~~~~~~~~~~~~~ ##
    #### Internal Functions ####
    ## ~~~~~~~~~~~~~~~~~~~~~~ ##
    
    ## Rescale the factor loadings when a Heywood case occurs
    commRescale <- function(lambda,
                            threshold = 1.0,
                            rescaleTo = .98) {
      
      ## Rescale the factor pattern matrix when Heywood case is detected
      # lambda is an unrotated factor structure matrix
      hsq.new <- hsq.old <- apply(lambda^2, 1, sum)
    
      ## communalities that exceed threshold will be rescaled
      hsq.new[which(hsq.old >= threshold)] <- rescaleTo
      
      ## Norm lambda by dividing each row by its communality
      ## Then pre-multiply by the new (sqrt) communality diagonal matrix
      newLambda <- diag( sqrt( hsq.new / hsq.old) ) %*% lambda
      
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
    
    ## If digits is not specified, give default values
    if ( is.null(digits) ) {
      digits <- options()$digits
    } # END if ( is.null(digits) )
    
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ##         First-level hierarchical factor analysis         ##
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    
    ## ----First Level ----
    
    ## ----____Extract Factors----
    
    ## Extract the factor solution
    L1Out <- faX(R          = R,
                 numFactors = numFactors[1],
                 facMethod  = facMethod,
                 faControl  = faControl)
    
    ## Save the factor loadings for rotation
    L1 <- L1Out$loadings
    
    ## Save level-one communality values
    HSquare <- L1Out$h2
    
    ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## Check for Heywood cases ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
    
    HeywoodFlag <- 0
    ## ----____Check for Heywood cases----
    
    ## If any communalities are greater than lambda > 1.0, re-scale problem items
    if ( L1Out$faFit$Heywood == TRUE ) {
      
    ## Rescale the communalities larger than 'threshold' to 'rescaleTo'
      L1 <- commRescale(lambda    = L1,
                        threshold = 1.0,
                        rescaleTo = rescaleH2)
      HeywoodFlag <- 1
     
      ## Find new communality values after rescaling
      HSquare <- apply(L1^2, 1, sum)
      
    } # END if ( L1Out$Heywood == TRUE )
    
    ## ----____Rotate first order solution----
    
    ## Rotate the factor structure
    firstRot <- faMain(urLoadings    = L1,
                       rotate        = rotate,
                       rotateControl = rotateControl)
    
    ## Store the Phi matrix and rotated loadings matrix
    L1 <- firstRot$loadings   ## Factor pattern
    R1 <- firstRot$Phi        ## Factor correlations
    
    ## Reorder the factors
    facOrder <- 
      orderFactors(Lambda  = L1,
                   PhiMat  = R1,
                   salient = .25,
                   reflect = TRUE)
    
    ## Reordered matrices
    L1 <- facOrder$Lambda
    R1 <- facOrder$PhiMat
    
    ## CG (11 Oct 2019): No need to retain uniquenesses, just need their sqrt
    ## Save sqrt of U_1 for later orthogonalization
    U_1 <- sqrt( diag(1 - as.vector(HSquare)) )
    
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ##         Second-level hierarchical factor analysis        ##
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    
    ## ----Second Level----
    
    ## ----____Extract Factors-----
    ## Unrotated factor loadings of the 2nd level of analysis
    L2Out <- faX(R          = R1,
                 numFactors = numFactors[2],
                 facMethod  = facMethod,
                 faControl  = faControl)
    
    ## Save the factor loadings
    L2 <- L2Out$loadings[]

    ## Extract communalities to create uniqueness matrix
    HSquare.2 <- L2Out$h2
    
    ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## Check for Heywood cases ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
    
    ## ----____Check for Heywood cases----
    
    ## If any communalities in L2 are greater than lambda > 1.0, re-scale problem items
    if ( L2Out$faFit$Heywood == TRUE ) {
      
      ## Rescale the communalities larger than 'threshold' to 'rescaleTo'
      L2 <- commRescale(lambda    = L2,
                        threshold = 1.0,
                        rescaleTo = rescaleH2)
      HeywoodFlag <- 2
      
      
      ## Find new communality values after rescaling
      HSquare.2.new <- apply(L2^2, 1, sum)
      
    } # END if ( L1Out$Heywood == TRUE )
    
    
    ## ----____Rotate or reflect second order solution----
    
    
    ## Only needs to rotate obliquely if more than 1 factor is extracted
    if (numFactors[2] > 1) {
      
      ## Rotate the factor structure
      secondRot <- faMain(urLoadings    = L2,
                          rotate        = rotate,
                          rotateControl = rotateControl)
      
      ## Store the Phi matrix and rotated loadings matrix
      L2 <- secondRot$loadings
      
      ## CG (11 Oct 2019): Used to retain uniquenesses, only need their sqrt
      ## Save for orthogonalization later
      U_2 <- sqrt( diag(1 - as.vector(HSquare.2)) )
      
      ## Save the phi matrix
      R2 <- secondRot$Phi
      
      ## Reorder factors
      facOrder2 <- 
        orderFactors(Lambda  = L2,
                     PhiMat  = R2,
                     salient = .25,
                     reflect = TRUE)
      
      ## Retain the ordered factor matrices
      L2 <- facOrder2$Lambda
      R2 <- facOrder2$PhiMat
      
    } else {  #END if(numFactors[2] > 1)
      
      ## CG (11 Oct 2019): Did not previously reflect in 1 gen fac models
      ## Reflect loadings of 1 factor model
      reflect <- sign( sum( L2 ) )
      
      ## CG (11 Oct 2019): Did not previously reflect in 1 gen fac models
      ## If reflection is 0, set to 1 
      if (reflect == 0) reflect <- 1
      
      ## CG (11 Oct 2019): Did not previously reflect in 1 gen fac models
      ## Orient the single factor
      L2 <- L2 * reflect
      
      ## CG (11 Oct 2019): Used to retain uniquenesses, only need their sqrt
      ## Save for orthogonalization later
      U_2 <- sqrt( diag(1 - as.vector(HSquare.2)) )
      
      ## No phi matrix when only 1 factor is present, set at unity
      R2 <- 1
      
    } # END if (numFactors[2] > 1)
    
    
    ## ---- 2-Level Orthogonalization ----
    B.out <- cbind(L1 %*% L2, L1 %*% U_2)
    
    ## Reflect all factors by default
    BSign <- diag( sign( colSums(B.out) ) )
    B.out <- B.out %*% BSign
    

    
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ##         Third-level hierarchical factor analysis         ##
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    
    ## If a 3rd level of hierarchy exists, begin the process
    
    ## ----Third Level-----
    
    ## ----____Extract Factors----
    
    if (length(numFactors) == 3) {
      
      ## Unrotated factor analysis of the 3rd order factor model
      L3Out <- faX(R          = R2,
                   numFactors = numFactors[3],
                   facMethod  = facMethod,
                   faControl  = faControl)
      
      ## Save the factor loadings
      L3 <- L3Out$loadings
      
      ## Compute the communality estimate
      HSquare.3 <- apply(L3^2, 1, sum)
      
      ## ----____Check for Heywood cases ----
      
      ## Detect if Heywood cases are present
      if ( L3Out$faFit$Heywood == TRUE ) {
        
        ## Rescale the communalities larger than 'threshold' to 'rescaleTo'
        L3 <- commRescale(lambda    = L3,
                          threshold = 1.0,
                          rescaleTo = rescaleH2)
        
        HeywoodFlag <- 3
        
        ## Find rescaled communality
        HSquare.3 <- apply(L3^2, 1, sum)
        
      } # END if ( any(HSquare.3 > 1.0) )
      
      ## ----____Rotate or reflect third order solution----
      
      ## If multiple factors exist in the 3rd level, rotate obliquely
      if (numFactors[3] > 1) {
        
        ## Conduct oblique rotation if multiple 3rd-order factors exist
        thirdRot <- faMain(urLoadings    = L3,
                           rotate        = rotate,
                           rotateControl = rotateControl)
        
        ## Obliquely rotated factor loadings in the 3rd level
        L3 <- thirdRot$loadings
        
        ## Set the 3rd level phi matrix
        R3 <- thirdRot$Phi
        
        ## Faster to compute from saved output
        U_3 <- sqrt( diag(1 - as.vector(HSquare.3)) )
        
      } else { ## if(numFactors[3] == 1)
        
        ## CG (11 Oct 2019): Did not previously reflect in 1 gen fac models
        ## Reflect loadings of 1 factor model
        reflect <- sign( sum( L3 ) )
        
        ## CG (11 Oct 2019): Did not previously reflect in 1 gen fac models
        ## If reflection is 0, set to 1 
        if (reflect == 0) reflect <- 1
        
        ## CG (11 Oct 2019): Did not previously reflect in 1 gen fac models
        ## Orient the single factor
        L3 <- L3 * reflect
        
        ## Compute the sqrt of the uniqueness matrix
        U_3 <- sqrt( diag(1 - as.vector(L3^2)) )
        
        ## Set phi matrix to unity (only 1 factor)
        R3 <- 1
        
      } # END if(numFactors[3] == 1)
      
      ## ---- 3-Level Orthogonalization ----
      
      ## Start orthogonalization for 3 levels
      B3    <- cbind(L3, U_3)
      B2    <- cbind((L2 %*% B3), U_2)
      B.out <- L1 %*% B2
      
      ## Reflect all factors by default
      BSign <- diag( sign( colSums(B.out) ) )
      B.out <- B.out %*% BSign
      
      ## ----____Dim names of 3 level solution ----
      
      ## Give the 3-level analysis some labels
      colnames(B.out) <-
        c(paste(rep("g", numFactors[3]), c(1:numFactors[3]), sep = ""),
          paste(rep("F", numFactors[2]), c(1:numFactors[2]), sep = ""),
          paste(rep("f", numFactors[1]), c(1:numFactors[1]), sep = ""))
      
      ## Names of items resemble that of initial correlation matrix <-
      rownames(B.out) <- rownames(R)
      # END if(length(numFactors) == 3))
      
    } else { ## Do this if only 2 levels exist
      
      ## ----____Dim names of 2 level solution ----
      
      ## Give the 2-level analysis some labels, ARBITRARY order
      colnames(B.out) <-
        c(paste(rep("g", numFactors[2]), c(1:numFactors[2]), sep = ""),
          paste(rep("F", numFactors[1]), c(1:numFactors[1]), sep = ""))
      
      ## Names of items resemble that of initial correlation matrix
      rownames(B.out) <- rownames(R)
      
    }## END if(length(numFactors) == 3){
    
    # Rescuale to maintain communalities of the L1 solution
     HSquare.B <- diag(B.out %*% t(B.out))
     B.out <- diag(sqrt(HSquare/HSquare.B)) %*% B.out
     
     if(HeywoodFlag > 0){
      cat("\n\nWARNING:  A Heywood case was detected i one of the factor solutions.\n")
     } 
     
    
    ## ----Available output----
    SLOutput <-
      list(L1            = round(L1, digits),    ## level 1 loadings
           L2            = round(L2, digits),    ## level 2 loadings
           L3            = NULL,  ## level 3 loadings
           Phi1          = round(R1, digits),    ## level 1 phi matrix
           Phi2          = round(R2, digits),    ## level 2 phi matrix
           Phi3          = NULL,  ## level 3 phi matrix
           U1            = round(U_1, digits),   ## Level 1 sqrt of uniqueness
           U2            = round(U_2, digits),   ## Level 2 sqrt of uniqueness
           U3            = NULL,  ## level 3 sqrt of uniquenesses
           B             = round(B.out, digits),
           rotateControl = rotateControl,
           faControl     = faControl,
           HeywoodFlag = HeywoodFlag) ## final bifactor solution
    
    ## If a 3-rd order solution is produced, populate the output of the 3rd level
    if ( length(numFactors) == 3 ) {
      
      SLOutput$L3   <- round(L3,  digits)
      SLOutput$Phi3 <- round(R3,  digits)
      SLOutput$Usq3 <- round(U_3, digits)
      
    } # END if ( length(numFactors) == 3 )
    
    ## Return the list of output
    SLOutput
    
  }## END SchmidLeiman()
