#' Generate Fungible Correlation Matrices
#' 
#' Generate fungible correlation matrices. For a given vector of standardized
#' regression coefficients, Beta, and a user-define R-squared value, Rsq, find
#' predictor correlation matrices, R, such that Beta' R Beta = Rsq. The size of
#' the smallest eigenvalue (Lp) of R can be defined.
#' 
#' 
#' @param R A p x p predictor correlation matrix.
#' @param Beta A p x 1 vector of standardized regression coefficients.
#' @param Lp Controls the size of the smallest eigenvalue of RstarLp.
#' @param eps Convergence criterion.
#' @param Print.Warnings Logical, default = TRUE. When TRUE, convergence
#' failures are printed.
#' @return \item{R}{Any input correlation matrix that satisfies Beta' R Beta =
#' Rsq.} \item{Beta}{Input vector of std reg coefficients.} \item{Rstar}{A
#' random fungible correlation matrix.} \item{RstarLp}{A fungible correlation
#' matrix with a fixed minimum eigenvalue (RstarLp can be PD, PSD, or ID).}
#' \item{s}{Scaling constant for Rstar.} \item{sLp}{Scaling constant for
#' RstarLp.} \item{Delta}{Vector in the null space of vecp(Beta Beta').}
#' \item{Q}{Left null space of Beta.} \item{FrobNorm}{Frobenius norm ||R -
#' Rstar||_F.} \item{FrobNormLp}{Frobenius norm ||R - RstarLp||_F given random
#' Delta.} \item{converged}{An integer code. 0 indicates successful
#' completion.}
#' @author Niels Waller
#' @references Waller, N. (2016). Fungible Correlation Matrices: A method for
#' generating nonsingular, singular, and improper correlation matrices for
#' Monte Carlo research. Multivariate Behavioral Research.
#' @keywords fungible
#' @export
#' @examples
#' 
#' 
#' library(fungible)
#' 
#' ## ===== Example 1 =====
#' ## Generate 5 random PD fungible R matrices
#' ## that are consistent with a user-defined predictive 
#' ## structure: B' Rxx B = .30
#' 
#' set.seed(246)
#' ## Create a 5 x 5 correlation matrix, R,  with all r_ij = .25
#' R.ex1 <- matrix(.25, 5, 5)
#' diag(R.ex1) <- 1
#' 
#' ## create a 5 x 1 vector of standardized regression coefficients, 
#' ## Beta.ex1 
#' Beta.ex1 <- c(-.4, -.2, 0, .2, .4)
#' cat("\nModel Rsq = ",  t(Beta.ex1) %*% R.ex1 %*% Beta.ex1)
#' 
#' ## Generate fungible correlation matrices, Rstar, with smallest
#' ## eigenvalues > 0.
#' 
#' Rstar.list <- list(rep(99,5)) 
#' i <- 0
#' while(i <= 5){
#'   out <- fungibleR(R = R.ex1, Beta = Beta.ex1, Lp = 1e-8, eps = 1e-8, 
#'                    Print.Warnings = TRUE)
#'   if(out$converged==0){
#'     i <- i + 1
#'     Rstar.list[[i]] <- out$Rstar
#'   }
#' }
#' 
#' ## Check Results
#' cat("\n *** Check Results ***")
#' for(i in 1:5){
#'   cat("\n\n\n+++++++++++++++++++++++++++++++++++++++++++++++++")
#'   cat("\nRstar", i,"\n")
#'   print(round(Rstar.list[[i]], 2),)
#'   cat("\neigenvalues of Rstar", i,"\n")
#'   print(eigen(Rstar.list[[i]])$values)
#'   cat("\nBeta' Rstar",i, "Beta = ",  
#'       t(Beta.ex1) %*% Rstar.list[[i]] %*% Beta.ex1)
#' }  
#' 
#' 
#' 
#' ## ===== Example 2 =====
#' ## Generate a PD fungible R matrix with a fixed smallest 
#' ## eigenvalue (Lp).
#' 
#' ## Create a 5 x 5 correlation matrix, R,  with all r_ij = .5
#' R <- matrix(.5, 5, 5)
#' diag(R) <- 1
#' 
#' ## create a 5 x 1 vector of standardized regression coefficients, Beta, 
#' ## such that Beta_i = .1 for all i 
#' Beta <- rep(.1, 5)
#' 
#' ## Generate fungible correlation matrices (a) Rstar and (b) RstarLp.
#' ## Set Lp = 0.12345678 so that the smallest eigenvalue (Lp) of RstarLp
#' ## = 0.12345678
#' out <- fungibleR(R, Beta, Lp = 0.12345678, eps = 1e-10, Print.Warnings = TRUE)
#' 
#' ## print R
#' cat("\nR: a user-specified seed matrix")
#' print(round(out$R,3)) 
#' 
#' ## Rstar
#' cat("\nRstar: A random fungible correlation matrix for R")
#' print(round(out$Rstar,3)) 
#' 
#' cat("\nCoefficient of determination when using R\n")
#' print(  t(Beta) %*% R %*% Beta )
#' 
#' cat("\nCoefficient of determination when using Rstar\n")
#' print( t(Beta) %*% out$Rstar %*% Beta)
#' 
#' ## Eigenvalues of  R
#' cat("\nEigenvalues of R\n")
#' print(round(eigen(out$R)$values, 9)) 
#' 
#' ## Eigenvalues of  Rstar
#' cat("\nEigenvalues of Rstar\n")
#' print(round(eigen(out$Rstar)$values, 9)) 
#' 
#' ## What is the Frobenius norm (Euclidean distance) between
#' ## R and Rstar
#' cat("\nFrobenious norm ||R - Rstar||\n")
#' print( out$FrobNorm)
#' 
#' ## RstarLp is a random fungible correlation matrix with 
#' ## a fixed smallest eigenvalue of 0.12345678
#' cat("\nRstarLp: a random fungible correlation matrix with a user-defined
#' smallest eigenvalue\n")
#' print(round(out$RstarLp, 3)) 
#' 
#' ## Eigenvalues of RstarLp
#' cat("\nEigenvalues of RstarLp")
#' print(eigen(out$RstarLp)$values, digits = 9) 
#' 
#' cat("\nCoefficient of determination when using RstarLp\n")
#' print( t(Beta) %*% out$RstarLp %*% Beta)
#' 
#' ## Check function convergence
#' if(out$converged) print("Falied to converge")
#' 
#' 
#' ## ===== Example 3 =====
#' ## This examples demonstrates how fungibleR  can be used
#' ## to generate improper correlation matrices (i.e., pseudo 
#' ## correlation matrices with negative eigenvalues).
#' library(fungible)
#' 
#' ## We desire an improper correlation matrix that
#' ## is close to a user-supplied seed matrix.  Create an 
#' ## interesting seed matrix that reflects a Big Five 
#' ## factor structure.
#' 
#' set.seed(123)
#' minCrossLoading <- -.2
#' maxCrossLoading <-  .2
#' F1 <- c(rep(.6,5),runif(20,minCrossLoading, maxCrossLoading))
#' F2 <- c(runif(5,minCrossLoading, maxCrossLoading), rep(.6,5), 
#'       runif(15,minCrossLoading, maxCrossLoading))
#' F3 <- c(runif(10,minCrossLoading,maxCrossLoading), rep(.6,5), 
#'       runif(10,minCrossLoading,maxCrossLoading) )
#' F4 <- c(runif(15,minCrossLoading,maxCrossLoading), rep(.6,5), 
#'       runif(5,minCrossLoading,maxCrossLoading))
#' F5 <- c(runif(20,minCrossLoading,maxCrossLoading), rep(.6,5))
#' FacMat <- cbind(F1,F2,F3,F4,F5)
#' R.bfi <- FacMat %*% t(FacMat)
#' diag(R.bfi) <- 1
#' 
#' ## Set Beta to a null vector to inform fungibleR that we are 
#' ## not interested in placing constraints on the predictive structure 
#' ## of the fungible R matrices. 
#' Beta <- rep(0, 25)
#' 
#' 
#' ## We seek a NPD fungible R matrix that is close to the bfi seed matrix.
#' ## To find a suitable matrix we generate a large number (e.g., 50000) 
#' ## fungible R matrices. For illustration purposes I will set Nmatrices
#' ## to a smaller number: 10.
#' Nmatrices<-10
#' 
#' ## Initialize a list to contain the Nmatrices fungible R objects
#' RstarLp.list <- as.list( rep(0, Nmatrices ) )
#' ## Initialize a vector for the Nmatrices Frobeius norms ||R - RstarLp||
#' FrobLp.vec <- rep(0, Nmatrices)
#' 
#' 
#' ## Constraint the smallest eigenvalue of RStarLp by setting
#' ## Lp = -.1 (or any suitably chosen user-defined value).
#' 
#' ## Generate Nmatrices fungibleR matrices and identify the NPD correlation 
#' ## matrix that is "closest" (has the smallest Frobenious norm) to the bfi 
#' ## seed matrix.
#' BestR.i <- 0
#' BestFrob <- 99
#' i <- 0
#' 
#' set.seed(1)
#' while(i < Nmatrices){
#'   out<-fungibleR(R = R.bfi, Beta, Lp = -.1, eps=1e-10) 
#'   ## retain solution if algorithm converged
#'   if(out$converged == 0)
#'   { 
#'     i<- i + 1
#'   ## print progress  
#'     cat("\nGenerating matrix ", i, " Current minimum ||R - RstarLp|| = ",BestFrob)
#'     tmp <- FrobLp.vec[i] <- out$FrobNormLp #Frobenious Norm ||R - RstarLp||
#'     RstarLp.list[[i]]<-out$RstarLp
#'     if( tmp < BestFrob )
#'     {
#'       BestR.i <- i     # matrix with lowest ||R - RstarLp||
#'       BestFrob <- tmp  # value of lowest ||R - RstarLp||
#'     }
#'   }
#' }
#' 
#' 
#' 
#' # CloseR is an improper correlation matrix that is close to the seed matrix. 
#' CloseR<-RstarLp.list[[BestR.i]]
#' 
#' plot(1:25, eigen(R.bfi)$values,
#'      type = "b", 
#'      lwd = 2,
#'      main = "Scree Plots for R and RstarLp",
#'      cex.main = 1.5,
#'      ylim = c(-.2,6),
#'      ylab = "Eigenvalues",
#'      xlab = "Dimensions")
#' points(1:25,eigen(CloseR)$values,
#'        type = "b",
#'        lty = 2,
#'        lwd = 2,
#'        col = "red")
#' 	   abline(h = 0, col = "grey")
#' legend(legend=c(expression(paste(lambda[i]~" of R",sep = "")),
#'                 expression(paste(lambda[i]~" of RstarLp",sep = ""))),
#'        lty=c(1,2),
#'        x = 17,y = 5.75,
#'        cex = 1.5,
#'        col=c("black","red"),
#'        text.width = 5.5,
#'        lwd = 2)
#' 
fungibleR<-function(R, Beta, Lp = .00, eps=1e-8, Print.Warnings=TRUE){
  #---------------------------------------------------------------#
  # fungibleR
  # Niels Waller
  # July 23, 2013
  # updated August 14, 2015
  # updated September 19, 2015
  # updated July 23, 2016
  #
  # Arguments:
  #       R            : a p x p predictor correlation matrix
  #       Beta         : a p x 1 vector of standardized regression 
  #                      coefficients
  #       Lp           : Controls the size of the
  #                      smallest eigenvalue of RstarLp
  #       eps          : convergence criterion
  #   Print.Warnings   : logical, default = TRUE. When TRUE, convergence failures are printed
  #
  # Returns:
  #       R            : input correlation matrix
  #       Beta         : input vector of std reg coefficients
  #       Rstar        : a random fungible correlation matrix
  #       RstarLp      : a fungible correlation matrix with a fixed minimum eigenvalue
  #                      (RstarLp can be PD, PSD, or NPD)
  #       s            : scaling constant for Rstar
  #       sLp          : scaling constant for RstarLp
  #       Delta        : vector in the null space of vecp(Beta Beta^T)
  #       FrobNorm     : Frobenius norm ||R - Rstar||_F 
  #       FrobNormLp   : Frobenius norm ||R - RstarLp||_F given random Delta
  #       converged    : An integer code. 0 indicates successful convergence. 
  #---------------------------------------------------------------# 

  p <- length(Beta)
  

  max.iterations <- 100
  
# Fnc: tr -- matrix trace
  tr <- function(X) sum(diag(X))
  
# Fnc: normF-- Frobenius norm
  normF <- function(X)   sqrt( tr(t(X) %*% X) )  

# Fnc: vecp
  vecp <- function(X){
    if(dim(X)[1]!=dim(X)[2]) stop( "matrix not square" )
    X[upper.tri( X, diag=FALSE )]       
  }

# Fnc: vecpInv -- inverse of vecp
  vecpInv <- function(x,p=2,X=NULL) {
    if(is.null(X)) X <- matrix(0,p,p)
    X[upper.tri(X,diag=FALSE)] <- x
    X + t(X)
  }

# vector of coefficient cross products
  b <- vecp(outer( Beta, Beta ))

# generate Z, a p x p-1 matrix of random deviates
  nNull <- (p * (p-1)/2) - 1  # dimension of null space
  
# When generating VERY large R matrices, restrict Delta 
# to a manageable subspace
  if(nNull > 25000){
    Z <- matrix(rnorm( (nNull+1)*100 ), nNull+1, 100)
    nNull <- 100
  }
  else
    Z <- matrix(rnorm( (nNull+1)*nNull ),nNull+1,nNull)

# append b to Z
  bZ<-cbind(b,Z)

# use QR decomposition to find null space of  b 
# Qstar = Q[,-1]  
  Q <- qr.Q(qr(bZ))
  Qstar <- Q[,-1]
  
 
# w = weights for null space basis vectors
  w <- rnorm(nNull, 0, 1)
  w <- w/sqrt(sum(w^2)) 

# generate Delta
  Delta <- Qstar %*% w
  Delta <- Delta/sqrt(sum(Delta^2)) # norm delta
 

# initialize s from Rayleigh Coefficient
  LpH <- eigen(vecpInv(Delta,p))$val[p]
  LpR <- eigen(R)$val[p]
  s <- s.init <-  (Lp - LpR)/LpH 
  

# for RstarLp fixed smallest eigenvalue  
# find sLp by Newton Raphson
# v.p is last eigenvector of (R + sH)
  H <- vecpInv(Delta,p=p)
  crit <- 99
  it <- 0
  while(crit >= eps){
    v.p<-eigen(R + s*H)$vectors[,p]
    s.old<-s
    s <- s.old -  as.numeric((eigen(R + s.old*H)$values[p] - Lp)/(2*(t(v.p)%*%H%*%v.p)))               
    crit<-abs(s-s.old) 
    it<-it+1
    if(it>max.iterations) {
             converged <- 1
             break()
    }
  }                  
    
  sLp<-s
  Rstar.Lp <- R + sLp*H
  

# for non fixed smallest eigenvalue
# find s for PSD Rstar 
# reinitialize s from Rayleigh Coefficient
  s <- s.init

  LpZero <- 0
# find s for LpZero = 0  by Newton Raphson
  crit <- 99
  it<-0
  while(crit>=eps){
    v.p<-eigen(R + s*H)$vectors[,p]
    s.old<-s
    s <- s.old - as.numeric((eigen(R + s.old*H)$values[p] - LpZero)/(2*(t(v.p)%*%H%*%v.p)))               
    crit<-abs(s-s.old) 
    it<-it+1
    if(it>max.iterations) {
      converged<-1
      break()
    }
  }                  
  
# bounds on s for random Rstar
  sLO<-s + 1e-8  # lower bound of s
  sHI <- (-LpR + 1)/LpH # upper bound of s
  sLOsHI <-sort( c(sLO, sHI) )
  
# generate a random PD Rstar 
  it<-0
  for(it in 1:max.iterations){
    s <- runif(1, sLOsHI[1], sLOsHI[2])
    Rstar <- R + s*H
    if(eigen(Rstar)$values[p] > 0) break()
  }  
  
  
# compute the Frobenius norm of the difference matrix  
  Fnorm <- normF(R - Rstar)  
  FnormLp <- normF(R - Rstar.Lp) 
  
  
  ##########################################################################
  ##            Convergence tests
  ##  1   Newton Raphson failed to converge in max.iterations
  ##  2   RstarLp correlation > 1 in absolute value
  ##  3   RstarLp minimum eigenvalue tolerance test failure
  ##  4   minimum eigenvalue of Rstar < 0   
  ##  5   Rstar correlation > 1 in absolute value
  ##  6   s is NA

# convergence test Rstar.Lp: Are all |r_ij|< 1 
  converged<-0 #initialize  
  Rcheck<-Rstar.Lp
  diag(Rcheck)<-0
  if( max(abs(Rcheck)) > 1) {
    if(Print.Warnings){
       warning(gettextf("Improper solution: Correlation > 1 in absolute value"), domain = NA)
    }  
    converged <- 2
  }
 # convergence test: Is minimum eigenvalue correct  
  if(abs(min(eigen(Rstar.Lp)$values) - Lp) > eps) {
    if(Print.Warnings){
       warning(gettextf("Minimum eigenvalue tolerance test failure"), domain = NA)
    }  
    converged <- 3
  }
  
  # convergence test: Is minimum eigenvalue of Rstar > 0   
  if(eigen(Rstar)$values[p] < 0)  converged <- 4
 
  
# convergence test Rstar: Are all |r_ij|< 1 
  Rcheck<-Rstar
  diag(Rcheck)<-0
  if( max(abs(Rcheck)) > 1) {
     if(Print.Warnings){
      warning(gettextf("Improper solution: Correlation > 1 in absolute value"), domain = NA)
    }  
    converged <- 5
  }
  
  if(is.na(s)) converged <- 6
  
  list( R = R,
        Beta = Beta,
        Rstar=Rstar, 
        RstarLp=Rstar.Lp,
        s = s,
        sLp=sLp,
        Delta=Delta, 
        Q = Q,
        FrobNorm=Fnorm,
        FrobNormLp=FnormLp,
        converged=converged)
} #END OF FUNCTION

#----------------------------------------------------------------


# Example 1


# set.seed(123)
# Beta <- c(.1,.2,.3,.4)
# R <- diag(4)
# Rsq <- t(Beta) %*% R %*% Beta
# Rsq
# out<-fungibleR(R,Beta,Lp = .12345, eps=1e-8)
# round(out$RstarLp,5)
# eigen(out$RstarLp)$values



# R <- matrix(.5,5,5)
# diag(R) <- 1
# Beta <- rep(.1,5)
# # 
# # #change value of Lp to control the size of the smallest eigenvalue of RstarMax
# out<-fungibleR(R,Beta,Lp = 0.12345678, eps=1e-8)
# cat("\nR\n")
# print(  round(out$R,3) )
# cat("\nRstar\n")
# print(  round(out$Rstar,3) )
# cat("\nEigen Rstar\n")
# print( round(eigen(out$Rstar)$values,9) )
# cat("\nRstarMax\n")
# print( round(out$RstarMax,3) )
# cat("\nEigen RstarMax\n")
# print(eigen(out$RstarMax)$values,digits=14) 
# print( t(Beta) %*% out$RstarMax %*% Beta )
# if(!out$converged) print("Falied to converge")


# Example 2 A large problm with >2000 predictors
#  set.seed(124)
#  Beta <- c(rep(.1,5), rep(0,2000), rep(.1,5))
#  R <- diag(2010)
#  Rsq <- t(Beta) %*% R %*% Beta
#  Rsq
#  out<-fungibleR(R,Beta,Lp = 0.01234500, eps=1e-12)
#  #print( round(out$RstarMax,3) )[1:100,1:100]
#  print( summary(out$RstarMax[lower.tri(R)]) )
#  print( eigen(round(out$RstarMax,7))$values)
# 
# 


