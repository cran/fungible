# Nov 4, 2017
# Construct single LS matrix
# This is 1000 times faster for 10 dimensions





#' Align the columns of two factor loading matrices
#' 
#' Align factor loading matrices across solutions using the Hungarian algorithm
#' to locate optimal matches. faAlign will match the factors of F2 (the input
#' matrix) to those in F1 (the target matrix) to minimize a least squares
#' discrepancy function or to maximize factor congruence coefficients (i.e.,
#' vector cosines).
#' 
#' 
#' @param F1 target Factor Loadings Matrix.
#' @param F2 input Factor Loadings Matrix. F2 will be aligned with the target
#' matrix, F1.
#' @param Phi2 optional factor correlation matrix for F2 (default = NULL).
#' @param MatchMethod "LS" (Least Squares) or "CC" (congruence coefficients).
#' @return \item{F2}{re-ordered and reflected loadings of F2.}
#' \item{Phi2}{reordered and reflected factor correlations.} \item{FactorMap}{a
#' 2 x k matrix (where k is the number of columns of F1) structured such that
#' row 1: the original column order of F2; row 2: the sorted column order of
#' F2.} \item{UniqueMatch}{(logical) indicates whether a unique match was
#' found.} \item{MatchMethod}{"LS" (least squares) or "CC" (congruence
#' coefficients, i.e., cosines).} \item{CC}{Congruence coefficients for the
#' matched factors.} \item{LS}{Root-mean-squared-deviations (least squares
#' criterion) for the matched factors.}
#' @note The Hungarian algorithm is implemented with the clue (Cluster
#' Ensembles, Hornik, 2005) package. See Hornik K (2005). A CLUE for CLUster
#' Ensembles. \emph{Journal of Statistical Software, 14}(12). doi:
#' 10.18637/jss.v014.i12 (URL: http://doi.org/10.18637/jss.v014.i12).
#' @author Niels Waller
#' @references Kuhn, H. W. (1955). The Hungarian Method for the assignment
#' problem. \emph{Naval Research Logistics Quarterly, 2}, 83-97.
#' 
#' Kuhn, H. W. (1956). Variants of the Hungarian method for assignment
#' problems. \emph{Naval Research Logistics Quarterly, 3}, 253-258.
#' 
#' Papadimitriou, C. & Steiglitz, K. (1982). Combinatorial Optimization:
#' Algorithms and Complexity. Englewood Cliffs: Prentice Hall.
#' @keywords fungible Statistics
#' @family Factor Analysis Routines
#' @import clue
#' @export
#' @examples
#' 
#' # This example demonstrates the computation of 
#' # non-parametric bootstrap confidence intervals
#' # for rotated factor loadings.
#' 
#' 
#' library(GPArotation)
#' 
#' data(HS9Var)
#' 
#' HS9 <- HS9Var[HS9Var$school == "Grant-White",7:15]
#' 
#' # Compute an R matrix for the HSVar9 Mental Abilities Data
#' R.HS9 <- cor(HS9)
#' 
#' varnames <- c( "vis.per", "cubes", 
#'             "lozenges", "paragraph.comp",
#'             "sentence.comp","word.mean",
#'             "speed.add", "speed.count.dots",
#'             "speed.discr")
#' 
#' 
#' 
#' # Extract and rotate a 3-factor solution
#' # via unweighted least squares factor extraction 
#' # and oblimin rotation. 
#' 
#' NFac <- 3
#' NVar <- 9
#' B <- 200      # Number of boostrap samples
#' NSubj <- nrow(HS9)
#' 
#' # Unrotated 3 factor uls solution 
#'  F3.uls <- fals(R = R.HS9, nfactors = NFac)
#'  
#' # Rotate via oblimin 
#'  F3.rot <- oblimin(F3.uls$loadings, 
#'                       gam = 0, 
#'                       normalize = FALSE)
#' 
#'  F3.loadings <- F3.rot$loadings
#'  F3.phi <- F3.rot$Phi
#'  
#'  # Reflect factors so that salient loadings are positive
#'  Dsgn <- diag(sign(colSums(F3.loadings^3)))
#'  F3.loadings <- F3.loadings %*% Dsgn
#'  F3.phi <- Dsgn %*% F3.phi %*% Dsgn
#'  
#'  rownames(F3.loadings) <- varnames
#'  colnames(F3.loadings) <- paste0("f", 1:3)
#'  colnames(F3.phi) <- rownames(F3.phi) <- paste0("f", 1:3)
#'  
#'  cat("\nOblimin rotated factor loadings for 9 Mental Abilities Variables")
#'  print( round(F3.loadings, 2))
#'  
#'  cat("\nFactor correlation matrix")
#'  print( round( F3.phi, 2))
#'  
#'   # Declare variables to hold bootstrap output
#'   Flist <- Philist <- as.list(rep(0, B))
#'   UniqueMatchVec <- rep(0, B)
#'   rows <- 1:NSubj
#'  
#'   # Analyze bootstrap samples and record results 
#'   for(i in 1:B){
#'     cat("\nWorking on sample ", i)
#'     set.seed(i)
#'     
#'     # Create bootstrap sanples
#'     bsRows <- sample(rows, NSubj, replace= TRUE)
#'     Fuls <- fals(R = cor(HS9[bsRows, ]), nfactors = NFac)
#'     # rotated loadings
#'     Fboot <- oblimin(Fuls$loadings,
#'                              gam = 0, 
#'                              normalize = FALSE)
#'     
#'     
#'     out <- faAlign(F1 = F3.loadings, 
#'                    F2 = Fboot$loadings, 
#'                    MatchMethod = "LS")
#'     
#'     Flist[[i]] <- out$F2 # aligned version of Fboot$loadings
#'     UniqueMatchVec[i] <- out$UniqueMatch
#'   }
#'   
#'   cat("\nNumber of Unique Matches: ", 
#'       100*round(mean(UniqueMatchVec),2),"%\n")
#' 
#'   
#'   #  Make a 3D array from list of matrices
#'   arr <- array( unlist(Flist) , c(NVar, NFac, B) )
#'   
#'   #  Get quantiles of factor elements over third dimension (samples)
#'   F95 <- apply( arr , 1:2 , quantile, .975 )
#'   F05 <- apply( arr , 1:2 , quantile, .025 )
#'   Fse <- apply( arr , 1:2, sd  )
#'   
#'   cat("\nUpper Bound 95% CI\n")
#'   print( round(F95,3))
#'   cat("\n\nLower Bound 95% CI\n")
#'   print( round(F05,3))
#'   
#'   # plot distribution of bootstrap estimates
#'   # for example element
#'   hist(arr[5,1,], xlim=c(.4,1),
#'        main = "Bootstrap Distribution for F[5,1]",
#'        xlab = "F[5,1]")
#'   
#'   print(round (F3.loadings, 2))
#'   cat("\nStandard Errors")
#'   print( round( Fse, 2))
#' 
faAlign <- function(F1, F2, Phi2 = NULL, MatchMethod = "LS"){
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ## date: Nov 4, 2017
  ## Author: Niels Waller
  ##
  ## faAlign will align the factors of F2 to F1
  ## by minimizing the least squares (LS) fit of (F1 - F2)
  ## or by maximizing the matched-factor congruence coefficients (CC).
  ##
  ## If a unique match is not found in the initial attempt then the program 
  ## uses the Hungarian algorithm to find an optimal match.
  ##
  ## see:
  ## Harold W. Kuhn, "The Hungarian Method for the assignment problem", Naval Research 
  ##   Logistics Quarterly, 2: 83-97, 1955. Kuhn's original publication.
  ## Harold W. Kuhn, "Variants of the Hungarian method for assignment problems", 
  ##   Naval Research Logistics Quarterly, 3: 253-258, 1956.
  ## C. Papadimitriou and K. Steiglitz (1982), Combinatorial Optimization:
  ##    Algorithms and Complexity. Englewood Cliffs: Prentice Hall.
  ##
  ## Requires the clue library
  ##
  ## Arguments:
  ##
  ## F1:      Order Factor Loadings Matrix
  ##
  ## F2:      Input Factor Loadings Matrix
  ##
  ## Phi2:    Optional factor correlation matrix for F2 (default = NULL)
  ##
  ##
  ## MatchMethod: "LS" (Least Squares) or "CC" (congruence coefficients)
  ##
  ##
  ## Values:
  ##
  ## F2 : The re-ordered and reflected loadings of F2
  ##
  ## Phi2: Reordered and reflected factor correlations
  ##
  ## FactorMap: a 2 x  Nfac matrix structured such that:
  ##    row 1: the original column order of F2
  ##    row 2: the sorted column order of F2  
  ##
  ## UniqueMatch: A logical indicating whether a unique match was found.
  ##
  ## MatchMethod: "LS" or "CC" (congruence coefficients)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  
 Nfac <- ncol(F1)

 UniqueMatch <- TRUE
 
  if(MatchMethod == "LS"){
  
     # Compute modified least squares (i.e., squared distance) matrix
     A <-  matrix(colSums(F1^2), Nfac, Nfac, byrow=FALSE)
     B <- t(matrix(colSums(F2^2),Nfac, Nfac, byrow=FALSE))
     
     # When factors are optimally reflected the cross product will be positive
     LSmat <- A + B  - 2 * abs(crossprod(F1,F2))
     
     # Test for unique matches
     Qmatch <- apply(LSmat,1,which.min)
     if(length(unique(Qmatch)) != Nfac){
       UniqueMatch <- FALSE
     }
   
     # If unique match not found minimize sum of squares  
     # library(clue)
     if(UniqueMatch == FALSE){  
       Qmatch <- clue::solve_LSAP(LSmat, maximum = FALSE)
     }
     
     FactorMap <- rbind(1:Nfac, Qmatch)
     
     # Find optimal reflections
     F1noZeros <- F1
     # This allows a column of F1 to have all zeros
     F1noZeros[,colSums(F1noZeros)==0]<-1
     
     Dsgn <- diag(sign(colSums( F1noZeros*F2[,Qmatch]) ))
     F2 <- F2[ , Qmatch] %*%  Dsgn

 }#End LS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## Match on Congruence Coefficients   
 if(MatchMethod == "CC"){
   D.F1 <- diag(1 / sqrt(apply(F1^2, 2, sum)))
   D.F2 <- diag(1 / sqrt(apply(F2^2, 2, sum)))

   # take absolute value of CC matrix
   absCosF1F2 <- abs( D.F1 %*% t(F1) %*% F2 %*% D.F2 )
   
   
   # Test for unique matches
   Qmatch <- apply(absCosF1F2, 1, which.max)
   if(length(unique(Qmatch)) != Nfac){
     UniqueMatch <- FALSE
   }
   
   # If unique match not found maximize avg CC  
   if(UniqueMatch == FALSE){  
      Qmatch <- clue::solve_LSAP(absCosF1F2, maximum = TRUE)
   }
   
   FactorMap <- rbind(1:Nfac, Qmatch)
   
   # Find optimal reflections
   Dsgn <- diag(sign(colSums( F1*F2[,Qmatch]) ))
   F2 <- F2[ , Qmatch] %*%  Dsgn
   
  
 }#End CC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

 rownames(FactorMap) <- c("Original Order", "Sorted Order")

 if(!is.null(Phi2)){
   Phi2 <- Dsgn %*% Phi2[Qmatch, Qmatch] %*% Dsgn
 }
 
 ## Compute Congrence Coefficients and 
 ## RMSE for final solution
 
 FIT <- function(F1,F2){
   D.F1 <- diag(1 / sqrt(apply(F1^2, 2, sum)))
   D.F2 <- diag(1 / sqrt(apply(F2^2, 2, sum)))
   
   # take absolute value of CC matrix
   CC <- diag( D.F1 %*% t(F1) %*% F2 %*% D.F2 )
   LS <- sqrt(apply( (F1 - F2)^2 , 2, mean) )
   list(CC = CC, LS = LS)
 }
 
 fit <- FIT(F1,F2)

 list(F2 = F2,
     Phi2 = Phi2,
     FactorMap = FactorMap,
     UniqueMatch = UniqueMatch,
     MatchMethod = MatchMethod,
     CC = fit$CC,
     LS = fit$LS)
} # End MatchFactors
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##




