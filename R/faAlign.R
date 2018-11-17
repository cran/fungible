# Nov 4, 2017
# Construct single LS matrix
# This is 1000 times faster for 10 dimensions



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
     Dsgn <- diag(sign(colSums( F1*F2[,Qmatch]) ))
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




