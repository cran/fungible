#' Factor Scores
#' 
#' This function computes factor scores by various methods. The function will acceptan an 
#' object of class  \code{faMain} or, alternatively,  user-input factor pattern (i.e., \code{Loadings}) and 
#' factor correlation (\code{Phi}) matrices.
#' 
#' @param X (Matrix) An N \code{x}  variables data matrix.  If X is a matrix of raw scores then
#'  \code{faScores} will convert the data to z scores. 
#' @param faMainObject (Object of class \strong{faMain}) The returned object 
#' from a call to \strong{faMain}. Default = NULL
#' @param Loadings (Matrix) A factor pattern matrix. Default = NULL.
#' @param Phi (Matrix) A factor correlation matrix.  Default = NULL. If a factor pattern 
#' is entered via the \strong{Loadings} argument but \strong{Phi = NULL} the program 
#' will set \code{Phi} to an identity matrix.
#' @param Method (Character) Factor scoring method. Defaults to the Thurstone or regression based method.
#' Available options include: 
#'   \itemize{
#'     \item{\strong{Thurstone}} Generates regression based factor score estimates.
#'     \item{\strong{Bartlett}}  Generates Bartlett method factor score estimates.
#'     \item{\strong{tenBerge}}  Generates factor score estimates with correlations identical 
#'     to that found in \strong{Phi}.
#'     \item{\strong{Anderson}} The Anderson Rubin method. Generates uncorrelated factor score estimates.  This 
#'     method is only appropriate for orthogonal factor models. 
#'     \item{\strong{Harman}} Generates estimated factor scores by Harman's idealized variables method.
#'     \item{\strong{PCA}} Returns unrotated principal component scores.
#'    }
#'  
#' @return 
#' \itemize{
#'     \item{\strong{fscores}} A matrix om common factor score estimates. 
#'     \item{\strong{Method}} The method used to create the factor score estimates.  
#'     \item{\strong{W}} The factor scoring coefficient matrix. 
#'     \item{\strong{Z}} A matrix of standardized data used to create the estimated factor scores.
#' }
#' @details \strong{faScores} can be used to calculate estimated factor scores by various methods.  
#' In general, to calculate score estimates,  users must input a data matrix \strong{X} and either (a) 
#' an object of class \strong{faMain} or (b) a factor loadings matrix, \strong{Loadings} and 
#' an optional (for oblique models) factor correlation matrix \strong{Phi}. The one exception to this rule 
#' concerns scores for the principal components model.  To calculate unrotated PCA scores (i.e., when 
#' \strong{Method = "PCA"}) users need only enter a data matrix, \strong{X}.
#' 
#' @author  Niels Waller
#' @references 
#'  \itemize{ 
#'    \item{Bartlett, M. S. (1937). The statistical conception of 
#'          mental factors.British Journal of Psychology, 28,97-104.}
#'    \item{Grice, J.  (2001).  Computing and evaluating factor scores.
#'         Psychological Methods, 6(4), 430-450.}
#'     \item{Harman, H. H.  (1976).  Modern factor analysis.  
#'      University of Chicago press.} 
#'    \item{McDonald, R. P. and Burr, E. J.  (1967).  A Comparison of Four Methods of 
#'    Constructing Factor Scores.  Psychometrika, 32, 381-401.}
#'    \item{Ten Berge, J. M. F., Krijnen, W. P., Wansbeek, T., and Shapiro, A.  
#'    (1999).  Some new results on correlation-preserving factor scores 
#'    prediction methods.  Linear Algebra and its Applications, 
#'    289(1-3), 311-318.}  
#'    \item{Tucker, L.  (1971).  Relations of factor score estimates to their use.
#'     Psychometrika, 36, 427-436.}
#'  }       
#' 
#' @family Factor Analysis Routines
#' 
#' @examples
#' lambda.Pop <- matrix(c(.41, .00, .00,
#'                        .45, .00, .00,
#'                        .53, .00, .00,
#'                        .00, .66, .00,
#'                        .00, .38, .00,
#'                        .00, .66, .00,
#'                        .00, .00, .68,
#'                        .00, .00, .56,
#'                        .00, .00, .55),
#'                        nrow = 9, ncol = 3, byrow = TRUE)
#'  NVar <- nrow(lambda.Pop)
#'  NFac <- 3
#'
#' ## Factor correlation matrix
#' Phi.Pop <- matrix(.50, nrow = 3, ncol = 3)
#' diag(Phi.Pop) <- 1
#'
#' #Model-implied correlation matrix
#' R <- lambda.Pop %*% Phi.Pop %*% t(lambda.Pop)
#' diag(R) <- 1
#'
#'
#' #Generate population data to perfectly reproduce pop R
#' Out <- simFA( Model = list(Model = "oblique"),
#'              Loadings = list(FacPattern = lambda.Pop),
#'              Phi = list(PhiType = "user",
#'                         UserPhi = Phi.Pop),
#'              FactorScores = list(FS = TRUE,
#'                                  CFSeed = 1,
#'                                  SFSeed = 2,
#'                                  EFSeed = 3,
#'                                  Population = TRUE,
#'                                  NFacScores = 100),
#'              Seed = 1)
#'
#' PopFactorScores <- Out$Scores$FactorScores
#' X <- PopObservedScores <- Out$Scores$ObservedScores
#'
#'
#'fout <- faMain(X             = X,
#'               numFactors    = 3,
#'               facMethod     = "fals",
#'               rotate        = "oblimin")
#'
#'
#' print( round(fout$loadings, 2) )
#' print( round(fout$Phi,2) )
#'
#' fload <- fout$loadings
#' Phi <- fout$Phi
#' 
#'   fsOut <- faScores(X = X, 
#'                     faMainObject = fout, 
#'                     Method = "Thurstone")
#'
#'   fscores <- fsOut$fscores
#'
#'   print( round(cor(fscores), 2 ))
#'   print(round(Phi,2)) 
#'
#'  CommonFS <- PopFactorScores[,1:NFac]
#'  SpecificFS <-PopFactorScores[ ,(NFac+1):(NFac+NVar)]
#'  ErrorFS <-   PopFactorScores[ , (NFac + NVar + 1):(NFac + 2*NVar) ]
#'
#' print( cor(fscores, CommonFS) )
#' @importFrom utils tail
#' @export

  faScores <- function(X = NULL, faMainObject = NULL,  Loadings = NULL, Phi = NULL, Method = "Thurstone"){

 
    
    # Compute Matrix Powers
     matrixPower <- function(M, power){
        ULU <- eigen(M)
        ULU$vectors %*% diag(ULU$values^power) %*% t(ULU$vectors)
     }
    
    ##--------------------------------------------#
    ## -- Error Checking ----
     
     if(!is.null(faMainObject)){
       fout <- faMainObject
       R <- fout$R
       fload <- faMainObject$loadings
       Phi <-   faMainObject$Phi
     }
     
     if(!is.null(Loadings)){
       fload <- Loadings
       if(is.null(Phi)) Phi <- diag(ncol(Loadings))
       R <- cor(X)
     }
     
     if(is.null(X)) stop("\n\nA data matrix must be specified\n\n")
     
     
      #  Check if X is in z score form
       Z <- X
       if( !all( apply(Z, 2, sd) == 1) ) Z <- scale(Z)
    
      # Check if faMainObject is of the appropriate class
       if(!is.null(faMainObject)){
          stopifnot(inherits(fout, "faMain"))
       }
    
  
    ## ---- Compute Factor Scores----
       
       #---- ___Thurstone ----
       if(Method=="Thurstone"){
         # called exact regression by Grice 2001a p. 71
         #Gorsuch p. 262
         
         # Calculate the eigenvalues of R
         eigsR <- eigen(R)$values
         
         # If last eigenvalue close to zero smooth matrix BY
         # to increase stability of inverse
         if(tail(eigsR, n=1) < 1E-6){
           R <- smoothBY(R, const=.98)$RBY
         }
         W <-  solve(R) %*% fload %*% Phi
         fscores <- Z %*% W
       } ##END Thurstone
    
       #---- ___Bartlett/ML estimates ----
    if(Method=="Bartlett"){
      Uinv2<-diag(1/diag(R - fload %*% Phi %*% t(fload)))
      #Gorsuch p. 264
      
      if(tail(eigen(t(fload) %*% Uinv2 %*%fload)$values, n=1) < 1E-8){
        stop("\n\n\n****FATAL ERROR****\nMatrix Computationally Singular\nUnable to compute Bartlett factor score estimates")
      }
      W <- Uinv2 %*% fload %*% solve(t(fload) %*% Uinv2 %*%fload)
      fscores <- Z %*% W
    }
       
     #---- ___ten Berge ----
       if(Method=="tenBerge"){
         #McDonald 1981
         # ten Berge et al 1999
         #Grice 2001b p. 433
         PhiSqrt <-matrixPower(M = Phi, power = .5) 
         L <- fload %*% PhiSqrt 
         RinvSqrt <- matrixPower(M = R, power = -.5)
         tLRinvL <- t(L) %*% solve(R) %*% L
         
         if(tail(eigen(tLRinvL)$values, n=1) < 1E-8){
           stop("\n\n\n****FATAL ERROR****\nMatrix Computationally Singular\nUnable to compute ten Berge factor score estimates")
         }
         
         
         MinvSqrt <- matrixPower( tLRinvL, power = -.5)
         C <-  RinvSqrt %*% L %*% MinvSqrt
         W <- RinvSqrt %*% C %*% PhiSqrt
         fscores <- Z %*% W
       }
       
       #---- ___Anderson Rubin ----
    if(Method=="Anderson"){
      #Anderson Rubin
      #Gorsuch p. 265
      Uinv2<-diag(1/diag(R - fload %*% Phi %*% t(fload)))
      M <- (t(fload) %*% Uinv2 %*% R %*% Uinv2 %*% fload)
      
      if(tail(eigen(M)$values, n=1) < 1E-8){
        stop("\n\n\n****FATAL ERROR****\nMatrix Computationally Singular\nUnable to compute Anderson Rubin  factor score estimates")
      }
      
      MinvSqrt <-  matrixPower(M = M, power = -.5)
      W <- Uinv2 %*% fload %*% MinvSqrt
      fscores <- Z %*% W
    }
    
    
    #---- ___Harman ----
    if(Method=="Harman"){
      # Harman
      #Grice 2001b Eq 9
      floadtfload <- fload %*% t(fload)
      
       if(tail(eigen(floadtfload)$values, n=1) < 1E-8){
          W <- MASS::ginv(floadtfload) %*% fload
       }else{
         W <- solve(floadtfload) %*% fload
       }   

      fscores <- Z %*% W
    }
       
      
    #----___PCA----
       #return unrotated principal component scores
       if(Method == "PCA"){
         W <- eigen(R)$vectors
         fscores <- Z %*% W
       }
       
       
    #---- RETURN LIST ----   
    list(fscores = fscores,
         Method = Method,
         W = W,
         data = Z)
  } ## END faScores
  
  ##############
