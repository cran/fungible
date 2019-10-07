#' Generate Thurstone's Box Data From length, width, and height box measurements
#'
#' Generate data for Thurstone's 20 variable and 26 variable Box Study From length, width, and height box measurements.  
#'
#' @param XYZ (Matrix) Length, width, and height measurements for N boxes.   The Amazon Box data 
#'  can be accessed by calling \code{data(AmxBoxes)}. The Thurstone Box data (20 hypothetical boxes) 
#'  can be accessed by calling \code{data(Thurstone20Boxes)}.
#' @param BoxStudy (Integer) If BoxStudy = 20 then data will be generated for 
#'    Thurstone's classic 20 variable box problem. If BoxStudy = 26 then data will 
#'    be generated for Thurstone's 26 variable box problem. Default: \code{BoxStudy = 20}.
#' @param Reliability (Scalar [0, 1] ) The common reliability value for each 
#'    measured variable. Default: Reliability = .75.
#' @param ModApproxErrVar (Scalar [0, 1] ) The proportion of reliable 
#'    variance (for each variable) that is due to all minor common factors. 
#'    Thus, if \code{x} (i.e., error free length) has variance \code{var(x)} and  
#'     \code{ModApproxErrVar = .10}, then \code{var(e.ma)/var(x + e.ma) = .10}.
#' @param  SampleSize (Integer) Specifies the number of boxes to be sampled from 
#'   the population.  If \code{SampleSize = NULL} then measurements will be 
#'   generated for the original input box sizes. 
#' @param  NMinorFac (Integer) The number of minor factors to use while 
#'   generating model approximation error. Default: \code{NMinorFac = 50.}
#' @param  epsTKL (Numeric [0, 1])  A parameter of the 
#'    Tucker, Koopman, and Linn (1969) algorithm that controls the spread of the influence of the minor factors.
#'    Default: \code{epsTKL = .20.} 
#' @param  Seed (Integer)  Starting seed for box sampling.     
#' @param SeedErrorFactors (Integer)  Starting seed for the error-factor scores. 
#' @param SeedMinorFactors (Integer)  Starting seed for the minor common-factor scores. 
#' @param PRINT (Logical) If PRINT = TRUE then the computed reliabilites will 
#'    be printed. Default: \code{PRINT = FALSE}.  Setting \code{PRINT} to TRUE can be useful 
#'    when \code{LB = TRUE}.
#' @param LB (lower bound; logical)  If LB = TRUE then minimum box measurements will be 
#'              set to LBVal (inches) if they
#'              fall below 0 after adding measurement error. If LB = FALSE then negative 
#'              attribute values will not be modified. This argument has no effect
#'              on data that include model approximation error. 
#' @param LBVal (Numeric) If \code{LB = TRUE} then values in \code{BoxDataE} will be bounded 
#'     from below at \code{LBVal}.  This can be used to avoid negative or very small box 
#'     measurements.
#' @param Constant (Numeric) Optional value to add to all box measurements. 
#' Default: \code{Constant = 0}.   
#'   
#' @details This function can be used with the Amazon boxes dataset (\code{data(AmzBoxes)}) or with any collection  
#' of user-supplied scores on three variables. The Amazon Boxes data were downloaded from the 
#' \code{BoxDimensions} website: (\url{https://www.boxdimensions.com/}). These data contain 
#' length (x),  width (y), and height (z) measurements for 98 Amazon shipping boxes.  In his 
#' classical monograph on Multiple Factor Analysis (Thurstone, 1947) Thurstone describes two data sets 
#' (one that he created from fictitious data and a second data set that he created from actual box measurements) 
#' that  were used to illustrate topics in factor analysis. The  first (fictitious) data set is 
#'  known as the Thurstone Box problem (see Kaiser and Horst, 1975).  To create his data for the Box problem, 
#'  Thurstone constructed 20 nonlinear combinations of fictitious length, width, and height measurements. 
#'  \strong{Box20} variables:
#' \enumerate{
#'      \item   x^2
#'      \item   y^2
#'      \item   z^2
#'      \item   xy
#'      \item   xz
#'      \item   yz
#'      \item   sqrt(x^2 + y^2)
#'      \item   sqrt(x^2 + z^2)
#'      \item   sqrt(y^2 + z^2)
#'      \item   2x + 2y
#'      \item   2x + 2z
#'      \item   2y + 2z 
#'      \item   log(x)
#'      \item   log(y)
#'      \item   log(z)
#'      \item   xyz
#'      \item   sqrt(x^2 + y^2 + z^2)
#'      \item   exp(x)
#'      \item   exp(y)
#'      \item   exp(z)
#'    }
#'    
#'  The second Thurstone Box problem contains measurements on the following 26 functions of length, width, and height. 
#'   \strong{Box26} variables:
#'   \enumerate{
#'   \item x
#'   \item y
#'   \item z
#'   \item xy 
#'   \item xz
#'   \item yz 
#'   \item x^2 * y
#'   \item x * y^2
#'   \item x^2 * z
#'   \item x * z^ 2
#'   \item y^2 * z
#'   \item y * z^2
#'   \item x/y
#'   \item y/x
#'   \item x/z
#'   \item  z/x
#'   \item  y/z
#'   \item  z/y
#'   \item 2x + 2y
#'   \item 2x + 2z
#'   \item 2y + 2z
#'   \item sqrt(x^2 + y^2)
#'   \item sqrt(x^2 + z^2)
#'   \item sqrt(y^2 + z^2)
#'   \item xyz
#'   \item sqrt(x^2 + y^2 + z^2)
#'   }
#' @details Note that when generating unreliable data (i.e., variables with 
#' reliability values less than 1) and/or data with model error, 
#' \strong{SampleSize} must be greater than \strong{NMinorFac}.
#' @return 
#' \itemize{
#'   \item \strong{XYZ} The length (x), width (y), and height (z) measurements for the sampled boxes. 
#'          If \code{SampleSize = NULL} then \code{XYZ} contains the x, y, z values for the 
#'          original 98 boxes.
#'   \item \strong{BoxData} Error free box measurements.  
#'   \item \strong{BoxDataE} Box data with added measurement error. 
#'   \item \strong{BoxDataEME}  Box data with added (reliable) model approximation and (unreliable) measurement error.
#'   \item \strong{Rel.E}  Classical reliabilities for the scores in \code{BoxDataE}.
#'   \item \strong{Rel.EME}  Classical reliabilities for the scores in \code{BoxDataEME}.
#'   \item \strong{NMinorFac}  Number of minor common factors used to generate \code{BoxDataEME}.
#'   \item \strong{epsTKL}  Minor factor spread parameter for the Tucker, Koopman, Linn algorithm.
#'   \item \strong{SeedErrorFactors}  Starting seed for the error-factor scores.
#'   \item \strong{SeedMinorFactors} Starting seed for the minor common-factor scores. 
#'   }
#' 
#'
#' @author
#' Niels G. Waller (nwaller@umn.edu)
#'
#'
#' @family Factor Analysis Routines
#' 
#' @references 
#' 
#' Cureton, E. E. & Mulaik, S. A. (1975). The weighted varimax rotation and the 
#' promax rotation. Psychometrika, 40(2), 183-195. 

#' Kaiser, H. F. and Horst, P.  (1975).  A score matrix for Thurstone's box problem.  
#' Multivariate Behavioral Research, 10(1), 17-26.  
#' 
#' Thurstone, L. L.  (1947).  Multiple Factor Analysis.  Chicago: 
#' University of Chicago Press. 
#' 
#' Tucker, L. R., Koopman, R. F., and Linn, R. L.  (1969).  Evaluation of factor 
#' analytic research procedures by means of simulated correlation matrices. 
#' \emph{Psychometrika, 34}(4), 421-459.  
#'
#' @examples
#' 
#'   data(AmzBoxes)
#'   BoxList <- GenerateBoxData (XYZ = AmzBoxes[,2:4],
#'                               BoxStudy = 20,  
#'                               Reliability = .75,
#'                               ModApproxErrVar = .10,
#'                               SampleSize = 300, 
#'                               NMinorFac = 50,
#'                               epsTKL = .20,
#'                               Seed = 1,
#'                               SeedErrorFactors = 1,
#'                               SeedMinorFactors = 2,
#'                               PRINT = FALSE,
#'                               LB = FALSE,
#'                               LBVal = 1,
#'                               Constant = 0)
#'                                 
#'    BoxData <- BoxList$BoxData
#'    
#'    RBoxes <- cor(BoxData)
#'    fout <- faMain(R = RBoxes,
#'                  numFactors = 3,
#'                  facMethod = "fals",
#'                  rotate = "geominQ",
#'                  rotateControl = list(numberStarts = 100,
#'                                       standardize = "CM")) 
#'                                       
#'   summary(fout)  
#' @export

GenerateBoxData <- function(XYZ,
                       BoxStudy = 20,      
                       Reliability = .75,
                       ModApproxErrVar = .10,
                       SampleSize = NULL, 
                       NMinorFac = 50,
                       epsTKL = .20,
                       Seed = 1,
                       SeedErrorFactors = 2,
                       SeedMinorFactors = 3,
                       PRINT = FALSE,
                       LB = FALSE,
                       LBVal = 1,
                       Constant = 0
                       ){

  if(Constant != 0){
    XYZ <- XYZ + Constant
  }  
  
  set.seed(Seed)
  
  ## ----Generate Sample Data ----
  if(!is.null(SampleSize)){
    rows <- sample(x = 1:nrow(XYZ), 
                   size = SampleSize, 
                   replace = TRUE, 
                   prob = NULL)
    XYZ <- XYZ[rows, ]
  }else{     #END Generate Sample Data
    SampleSize <- nrow(XYZ)
  }  
  
  #----Error Checking----
  if(ModApproxErrVar > 0){
      if(SampleSize < NMinorFac){ 
        stop("\n*** FATAL ERROR: Insufficient sample size for model parameters ***")
      }
  }
  

  # x = length, 
  # y = width 
  # z = height
  # xT = true (input) x values 
  x <- xT <- XYZ[, 1]
  y <- yT <- XYZ[, 2]
  z <- zT <- XYZ[, 3]
  
  NBoxes <- length(x)
  if(BoxStudy == 20) NVar = 20 
  if(BoxStudy == 26) NVar = 26 
  
  
  
  ## ---- AddRandomError ----
  AddRandomError <- function(X, Rel, iter) {
    
    # X is a single variable
    # REL is the desired reliability
    # When X comes in it is error free
    V_T <- var(X)
    
    e.sd <- sqrt( (V_T / Rel) - V_T)
    
    random.error <- rnorm(NBoxes)
    e <- scale(resid(lm(random.error ~ X))) * e.sd
 
    # Xobs equal true plus error scores
    Xobs <- X + e
    
    # Dont allow negative measurements
    # minimum box dimension = 1 inch
    if ( isTRUE(LB) ) Xobs[Xobs <= 0] <- LBVal
    if(PRINT) cat("Var ", iter, " reliability  = ", round(V_T/var(Xobs),3),"\n")
    Xobs
  }#END  AddRandomError
  
  
  
  ##---- Minor Common Factors: Tucker Koopman Linn ----
    AddModelError <- function(X,
                            MEVar,
                            NMinorFac = 50,
                            epsTKL = .20){
    
    # MEVar comes in as a percent of total variance
    # Need to scale it
    gamma <- MEVar/(1 - MEVar)
    
    # X is input data (98 by 26)
    # MEVar = model approximation error
    
    # scale to simplify metric of model approximation error
    X <- scale(X)

    
    # Create matrix of minor factors
    # MacCallum and Tucker 1991 chose
    # 50 minor factors.
    # Briggs and MacCallum chose 150 minor factors
    NMinorFactors <- NMinorFac
    W <- matrix(0, ncol(X), NMinorFactors)
    
    # s is a scaling constant
    s <- 1 - epsTKL
    
    #generate unscaled random loadings
    for (i in 0:(NMinorFactors - 1)) {
      W[, i+1] <- rnorm(NVar, 
                        mean = 0, 
                        sd   = s^i)  #TKL EQ 13, p. 431
    }# END for (i in 0:(NMinorFactors - 1)) {
    
   
    # row SS of unscaled random loadings
    wsq <- diag(W %*% t(W))
    
    #
    ModelErrorVar <- rep(gamma, NVar)
  
    # Rescale loadings
    W <- diag(sqrt(ModelErrorVar/wsq)) %*% W
     

    # Create scores for minor factors
     ranData <- matrix(rnorm(NBoxes * NMinorFactors), 
                       nrow = NBoxes,
                       ncol = NMinorFactors)
     
     # these should be orthogonal to X
     for(i in 1:NMinorFactors){
       ranData[,i] <- scale(resid(lm(ranData[,i] ~ X)))
     }   
     
     # Orthogonalize and scale minor factor scores
     Z.ME <- scale( svd(ranData)$u ) 
 

     ModelErrorFactorScores <- Z.ME %*% t(W)
    
    # scale final scores to Z-scores
     BoxDataME <- scale(X + ModelErrorFactorScores)
    
    # check
    #print( var(X[,19])/var(X[,19] + ModelErrorFactorScores[,19]) ) 

    BoxDataME
    
  }#END AddModelError
  
  
  ##----Create Boxes ----
  # ---- Create attributes for Thurstone's 20 Boxes Example----
  MakeBox20 <- function(x, y, z){
    xsq <- x^2
    ysq <- y^2
    zsq <- z^2
    xy <- x * y
    xz <- x * z
    yz <- y * z
    root.xsq.plus.ysq <- sqrt(x^2 + y^2)
    root.xsq.plus.zsq <- sqrt(x^2 + z^2)
    root.ysq.plus.zsq <- sqrt(y^2 + z^2)
    twox.plus.twoy <- 2 * x + 2 * y
    twox.plus.twoz <- 2 * x + 2 * z
    twoy.plus.twoz <- 2 * y + 2 * z 
    logx <- log(x)
    logy <- log(y)
    logz <- log(z)
    xyz <- x * y * z
    root.xsq.plus.ysq.plus.zsq <- sqrt(x^2 + y^2 + z^2)
    expx <- exp(x)
    expy <- exp(y)
    expz <- exp(z)
    
    BoxData <- cbind(
      xsq,
      ysq,
      zsq,
      xy,
      xz,
      yz,
      root.xsq.plus.ysq,
      root.xsq.plus.zsq,
      root.ysq.plus.zsq,
      twox.plus.twoy,
      twox.plus.twoz,
      twoy.plus.twoz,
      logx,
      logy,
      logz,
      xyz,
      root.xsq.plus.ysq.plus.zsq,
      expx,
      expy,
      expz)
    
    BoxData
  }
  
  # ---- Create attributes for Thurstone's 26 Boxes Example----
  MakeBox26 <- function(x, y, z) {
    xy <- x * y
    xz <- x * z
    yz <- y * z
    xsqy <- x^2 * y
    xysq <- x * y^2
    xsqz <- x^2 * z
    xzsq <- x * z^2
    ysqz <- y^2 * z
    yzsq <- y * z^2
    x.over.y <- x/y
    y.over.x <- y/x
    x.over.z <- x/z
    z.over.x <- z/x
    y.over.z <- y/z
    z.over.y <- z/y
    twox.plus.twoy <- 2 * x + 2 * y
    twox.plus.twoz <- 2 * x + 2 * z
    twoy.plus.twoz <- 2 * y + 2 * z
    root.xsq.plus.ysq <- sqrt(x^2 + y^2)
    root.xsq.plus.zsq <- sqrt(x^2 + z^2)
    root.ysq.plus.zsq <- sqrt(y^2 + z^2)
    xyz <- x * y * z
    root.xsq.plus.ysq.plus.zsq <- sqrt(x^2 + y^2 + z^2)
    
    BoxData <- cbind(
      x,
      y,
      z,
      xy,
      xz,
      yz,
      xsqy,
      xysq,
      xsqz,
      xzsq,
      ysqz,
      yzsq,
      x.over.y,
      y.over.x,
      x.over.z,
      z.over.x,
      y.over.z,
      z.over.y,
      twox.plus.twoy,
      twox.plus.twoz,
      twoy.plus.twoz,
      root.xsq.plus.ysq,
      root.xsq.plus.zsq,
      root.ysq.plus.zsq,
      xyz,
      root.xsq.plus.ysq.plus.zsq
    )
    
    return(BoxData)
  }# END MakeIndicators
  
 
  # Error Free Box Data
    if(BoxStudy == 20){
       BoxData <- MakeBox20(xT, yT, zT)
    }
   if(BoxStudy == 26){
      BoxData <- MakeBox26(xT, yT, zT)
   }
  
   BoxDataE <- NULL
   # ----Create BoxDataE----
   if(Reliability < 1){
       BoxDataE <- BoxData
       set.seed(SeedErrorFactors)
       for (j in 1:ncol(BoxDataE)) {
         BoxDataE[, j] <- AddRandomError(BoxDataE[, j], 
                                           Reliability, 
                                           iter = j)
       }
   } # END  if(Reliability < 1) 
   
   
 
  # # ----Create BoxDataEME----
   BoxDataEME <- NULL
   if(ModApproxErrVar > 0){
       # First add model Error  
       set.seed(SeedMinorFactors)
       BoxDataME <- AddModelError(BoxData,
                                  MEVar = ModApproxErrVar,
                                  NMinorFac = NMinorFac,
                                  epsTKL = epsTKL)
       # Then add random error
       set.seed(SeedErrorFactors)
       BoxDataEME <- BoxDataME
       for (j in 1:ncol(BoxDataEME)) {
       BoxDataEME[, j] <- AddRandomError(BoxDataME[, j], 
                                         Reliability, 
                                         iter = j)
       }
   }# END if(ModApproxErrVar > 0)   
  
  # Reliability of Box data with error
  Rel.E <- Rel.EME <- rep(1, ncol(BoxData))
  if(Reliability < 1){
  Rel.E <- apply(BoxData, 2, var)/apply(BoxDataE, 2, var)
  }
  
  # Reliability of Box data with model approximation error and
  # random error
  Rel.EME <- NULL
  if(ModApproxErrVar > 0) {
      Rel.EME <- 1/apply(BoxDataEME, 2, var)
  }    
  
 
  list(XYZ = XYZ,
       BoxData = BoxData,
       BoxDataE = BoxDataE,
       BoxDataEME = BoxDataEME,
       Rel.E = Rel.E,
       Rel.EME = Rel.EME,
       NMinorFac = NMinorFac,
       epsTKL = epsTKL,
       SeedErrorFactors = SeedErrorFactors,
       SeedMinorFactors = SeedMinorFactors)
} #END CreateBoxData






