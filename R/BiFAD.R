#' Bifactor Analysis via Direct Schmid-Leiman (DSL) Transformations
#'
#' This function estimates the (rank-deficient) Direct Schmid-Leiman (DSL) bifactor solution as well as the (full-rank) Direct Bifactor (DBF) solution.
#'
#' @param R (Matrix) A correlation matrix.
#' @param B (Matrix) Bifactor target matrix. If B is NULL the program will create an empirically defined target matrix.
#' @param numFactors (Numeric) The number of group factors to estimate.
#' @param rotate (Character) Designate which rotation algorithm to apply. See the \code{\link{faMain}} function for more details about possible rotations. An oblimin rotation is the default.
#' @param salient (Numeric) Threshold value for creating an empirical target matrix.
#' @inheritParams faMain
#'
#'
#' @return The following output are returned in addition to the estimated Direct Schmid-Leiman bifactor solution.
#'
#' \itemize{
#'   \item \strong{B}: (Matrix) The target matrix used for the Procrustes rotation.
#'   \item \strong{BstarSL}: (Matrix) The resulting (rank-deficient) matrix of Direct Schmid-Leiman factor loadings.
#'   \item \strong{BstarFR}: (Matrix) The resulting (full-rank) matrix of Direct Bifactor factor loadings.
#'   \item \strong{rmsrSL}: (Scalar) The root mean squared residual (rmsr) between the known B matrix and the estimated (rank-deficient) Direct Schmid-Leiman rotation. If the B target matrix is empirically generated, this value is NULL.
#'   \item \strong{rmsrFR}: (Scalar) The root mean squared residual (rmsr) between the known B matrix and the estimated (full-rank) Direct Bifactor rotation. If the B target matrix is empirically generated, this value is NULL.
#' }
#'
#' @author
#' \itemize{
#'   \item Niels G. Waller (nwaller@umn.edu)
#' }
#'
#' @family Factor Analysis Routines
#' 
#' @references
#' \itemize{
#'   \item Giordano, C. & Waller, N. G. (under review). Recovering bifactor models: A comparison of seven methods.
#'   \item Mansolf, M., & Reise, S. P. (2016). Exploratory bifactor analysis: The Schmid-Leiman orthogonalization and Jennrich-Bentler analytic rotations. \emph{Multivariate Behavioral Research, 51}(5), 698-717.
#'   \item Waller, N. G. (2018). Direct Schmid Leiman transformations and rank deficient loadings matrices. \emph{Psychometrika, 83}, 858-870.
#' }
#'
#'
#' @examples
#' cat("\nExample 1:\nEmpirical Target Matrix:\n")
#' # Mansolf and Reise Table 2 Example
#' Btrue <- matrix(c(.48, .40,  0,   0,   0,
#'                   .51, .35,  0,   0,   0,
#'                   .67, .62,  0,   0,   0,
#'                   .34, .55,  0,   0,   0,
#'                   .44,  0, .45,   0,   0,
#'                   .40,  0, .48,   0,   0,
#'                   .32,  0, .70,   0,   0,
#'                   .45,  0, .54,   0,   0,
#'                   .55,  0,   0, .43,   0,
#'                   .33,  0,   0, .33,   0,
#'                   .52,  0,   0, .51,   0,
#'                   .35,  0,   0, .69,   0,
#'                   .32,  0,   0,   0, .65,
#'                   .66,  0,   0,   0, .51,
#'                   .68,  0,   0,   0, .39,
#'                   .32,  0,   0,   0, .56), 16, 5, byrow=TRUE)
#'
#' Rex1 <- Btrue %*% t(Btrue)
#' diag(Rex1) <- 1
#'
#' out.ex1 <- BiFAD(R          = Rex1,
#'                  B          = NULL,
#'                  numFactors = 4,
#'                  facMethod  = "fals",
#'                  rotate     = "oblimin",
#'                  salient    = .25)
#'
#' cat("\nRank Deficient Bifactor Solution:\n")
#' print( round(out.ex1$BstarSL, 2) )
#'
#' cat("\nFull Rank Bifactor Solution:\n")
#' print( round(out.ex1$BstarFR, 2) )
#'
#' cat("\nExample 2:\nUser Defined Target Matrix:\n")
#'
#' Bpattern <- matrix(c( 1,  1,  0,   0,   0,
#'                       1,  1,  0,   0,   0,
#'                       1,  1,  0,   0,   0,
#'                       1,  1,  0,   0,   0,
#'                       1,  0,  1,   0,   0,
#'                       1,  0,  1,   0,   0,
#'                       1,  0,  1,   0,   0,
#'                       1,  0,  1,   0,   0,
#'                       1,  0,   0,  1,   0,
#'                       1,  0,   0,  1,   0,
#'                       1,  0,   0,  1,   0,
#'                       1,  0,   0,  1,   0,
#'                       1,  0,   0,   0,  1,
#'                       1,  0,   0,   0,  1,
#'                       1,  0,   0,   0,  1,
#'                       1,  0,   0,   0,  1), 16, 5, byrow=TRUE)
#'
#' out.ex2 <- BiFAD(R          = Rex1,
#'                  B          = Bpattern,
#'                  numFactors = NULL,
#'                  facMethod  = "fals",
#'                  rotate     = "oblimin",
#'                  salient    = .25)
#'
#' cat("\nRank Deficient Bifactor Solution:\n")
#' print( round(out.ex2$BstarSL, 2) )
#'
#' cat("\nFull Rank Bifactor Solution:\n")
#' print( round(out.ex2$BstarFR, 2) )
#'
#' @export

BiFAD <- function(R,
                  B             = NULL,
                  numFactors    = NULL,
                  facMethod     = "fals",
                  rotate        = "oblimin",
                  salient       = .25,
                  digits        = NULL,
                  rotateControl = NULL,
                  faControl     = NULL) {
  
  ###############################################################
  ## AUTHOR: Niels Waller
  ## August 18, 2017
  ##  requires  psych:   package for initial factor extraction
  ##            GPArotation:  for rotation options
  ##
  ## Arguments:
  ##  R:  Input correlation matrix
  ##
  ##  B:  bifactor target matrix. If B=NULL the program will
  ##        create an empirically defined target matrix.
  ##
  ##  nGroup: Number of group factors in bifactor solution.
  ##
  ##  factorMethod: factor extraction method. Options include:
  ##        minres (minimum residual), ml (maximum likelihood),
  ##        pa (principal axis), gls (generalized least squares).
  ##
  ##  rotation:  factor rotation method. Current options include:
  ##        oblimin, geominQ, quartimin, promax.
  ##
  ##  salient:  Threshold value for creating an empirical target
  ##        matrix.
  ##
  ##  maxitFA:  Maximum iterations for the factor extraction
  ##        method.
  ##
  ##  maxitRotate: Maximum iterations for the gradient pursuit
  ##        rotation algorithm.
  ##
  ##  gamma: Optional tuning parameter for oblimin rotation.
  ##
  ##  Value:
  ##
  ##    B: User defined or empirically generated target matrix.
  ##
  ##    BstarSL:  Direct S-L solution.
  ##
  ##    BstarFR:  Direct full rank bifactor solution.
  ##
  ##    rmsrSL:   Root mean squared residual of (B - BstarSL).
  ##
  ##    rmsrFR:   Root mean squared residual of (B - BstarFR).
  ################################################################
  
  ## ~~~~~~~~~~~~~~~~~~ ##
  #### Error Checking ####
  ## ~~~~~~~~~~~~~~~~~~ ##
  
  ## FATAL ERROR
  if ( is.null(B) & is.null(numFactors) ) {
    stop("\n\n FATAL ERROR: Either B or numFactors must be specified")
  } # END if (is.null(B) & is.null(nGroup))
  
  ## Set default for digits if null
  if ( is.null(digits) ) {
    digits <- options()$digits
  } # END if ( is.null(digits) )
  
  ## ~~~~~~~~~~~~~~~~~~ ##
  #### Begin Function ####
  ## ~~~~~~~~~~~~~~~~~~ ##
  
  ## Compute orthogonal Procrustes rotation
  Procrustes <-function(M1, M2){
    tM1M2 <- t(M1) %*% M2
    svdtM1M2 <- svd(tM1M2)
    P <- svdtM1M2$u
    Q <- svdtM1M2$v
    T <- Q %*% t(P)
    ## Orthogonally rotate M2 to M1
    M2 %*% T
  } # END Procrustes
  
  ## Compute unrotated factor solution
  if ( is.null(numFactors) ) {
    numFactors <- ncol(B) - 1
  } # END if ( is.null(numFactors) )
  
  ## Extract rank-deficient factor loadings
  faLoad <- faX(R          = R,
                numFactors = numFactors,
                facMethod  = facMethod,
                faControl  = faControl)$loadings[]
  
  ## Append column of zeros to create rank deficient matrix
  L0 <- cbind(faLoad, 0)
  
  ## Create 0/1 Target matrix if B not supplied
  Bflag <- 1 # set to 1 if Target matrix (B) is known
  
  ## Create 0/1 Target matrix if B not supplied
  if ( is.null(B) )  Bflag <- 0
  
  if ( !Bflag & ( !is.null(numFactors) ) ) {
    
    ## Rotate the rank-deficient factor structure to find B target
    B <- faMain(urLoadings    = faLoad,
                rotate        = rotate,
                rotateControl = rotateControl)$loadings[]
    
    ## Record signs of loadings
    signB <- sign(B)
    
    ## Convert Target matrix into signed 0/1 matrix
    B <- signB * matrix(as.numeric(abs(B) >= salient),
                        nrow = nrow(R),
                        ncol = numFactors)
    
    ## Append ones vector for general factor
    B <- cbind(1, B)
    
  } ## END if ( !Bflag & ( !is.null(numFactors) ) ) {
  
  ## Rotate rank deficient loading matrix to Target
  BstarSL <- Procrustes(B, L0)
  
  ## Compute non-hierarchical bifactor solution
  F2 <- faX(R          = R,
            numFactors = numFactors + 1,
            facMethod  = facMethod,
            faControl  = faControl)$loadings[]
  
  ## Rotate full rank loading matrix to Target
  BstarFR <- Procrustes(B, F2)
  
  ## Compute root mean squared residual of Target matrix (B)
  ## and best fitting SL and FR solutions (Bstar)
  rmsrSL <- rmsrFR <- NA
  if (Bflag == 1) {
    rmsrSL = sqrt(mean((B - BstarSL)^2))
    rmsrFR = sqrt(mean((B - BstarFR)^2))
  } # END if (Bflag == 1)
  
  list(B       = round(B,       digits),
       BstarSL = round(BstarSL, digits),
       BstarFR = round(BstarFR, digits),
       rmsrSL  = rmsrSL,
       rmsrFR  = rmsrFR)
  
} ## END BiFAD

