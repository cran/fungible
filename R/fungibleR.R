
fungibleR<-function(R, Beta, Lp = .00, eps=1e-16){
  #---------------------------------------------------------------#
  # fungibleR
  # Niels Waller
  # July 23, 2013
  # updated August 14, 2015
  #
  # Arguments:
  #       R            : a p x p predictor correlation matrix
  #       Beta         : a p x 1 vector of standardized regression 
  #                      coefficients
  #       Lp           : Controls the size of the
  #                      smallest eigenvalue of RstarLp
  #       eps          : convergence criterion
  #
  # Returns:
  #       R            : input correlation matrix
  #       Beta         : input vector of std reg coefficients
  #       Rstar        : a random fungible correlation matrix
  #       RstarLp      : a fungible correlation matrix with a fixed minimum eigenvalue
  #                      (RstarLp can be PD, PSD, or ID)
  #       s            : scaling constant for Rstar
  #       sLp          : scaling constant for RstarLp
  #       Delta        : vector in the null space of vecp(Beta Beta^T)
  #       FrobNorm     : Frobenius norm ||R - Rstar||_F 
  #       FrobNormLp   : Frobenius norm ||R - RstarLp||_F given random Delta
  #       converged    : a logical indicating whether the algorithm has converged 
  #---------------------------------------------------------------# 

  p <- length(Beta)
  
  
# Fnc: tr -- matrix trace
  tr <- function(X) sum(diag(X))
  
# Fnc: normF-- Frobenius norm
  normF <- function(X)   sqrt(tr(t(X) %*% X))  

# Fnc: vecp
  vecp <- function(X){
    if(dim(X)[1]!=dim(X)[2]) stop("matrix not square")
    X[upper.tri(X,diag=FALSE)]       
  }

# Fnc: vecpInv -- inverse of vecp
  vecpInv <- function(x,p=2,X=NULL) {
    if(is.null(X)) X<-matrix(0,p,p)
    X[upper.tri(X,diag=FALSE)] <- x
    X + t(X)
  }

# vector of coefficient cross products
  b <- vecp(outer(Beta,Beta))

# generate Z, a p x p-1 matrix of random deviates
  nNull <- (p * (p-1)/2) - 1  # dimension of null space
  
# When generating VERY large R matrices, restrict Delta 
# to a manageable subspace
  if(nNull>25000){
    Z <- matrix(rnorm( (nNull+1)*100 ),nNull+1,100)
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
  w<- rnorm(nNull,0,1)
  w <- w/sqrt(sum(w^2))

# generate Delta
  Delta <- Qstar %*% w
  Delta <- Delta/sqrt(sum(Delta^2)) # norm delta
 

# initialize s from Rayleigh Coefficient
  LpH <- eigen(vecpInv(Delta,p))$val[p]
  LpR <- eigen(R)$val[p]
  s <- (Lp - LpR)/LpH 


# find sLp by Newton Raphson
# v.p is last eigenvector of (R + sH)
  H <- vecpInv(Delta,p=p)
  crit <- 99
  it<-0
  while(crit>=eps){
    v.p<-eigen(R + s*H)$vectors[,p]
    s.old<-s
    s <- s.old - as.numeric((eigen(R + s.old*H)$values[p] - Lp)/(2*(t(v.p)%*%H%*%v.p)))               
    crit<-abs(s-s.old) 
    it<-it+1
    if(it>500) {
             converged<-0
             break()
    }
  }                  
    
  sLp<-s
  Rstar.Lp <- R + sLp*H
  

# generate a random Rstar by multiplying 
# Delta by a random number between eps (a small number > 0) and sLp
  eps<-1e-6
  s<-runif(1,eps,sLp)
  Rstar <- R + s*H 

# compute the Frobenius norm of the difference matrix  
  Fnorm <- normF(R - Rstar)  
  FnormLp <- normF(R - Rstar.Lp) 

# convergence test: Are all |r_ij|< 1 
  converged<-0 #initialize  
  Rcheck<-Rstar.Lp
  diag(Rcheck)<-0
  if( max(abs(Rcheck)) > 1) {
    warning(gettextf("Improper solution: Correlation > 1 in absolute value"), domain = NA)
    converged <- 1
  }
  
  if(is.na(s)) converged<-1
  
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


#set.seed(123)
#Beta <- c(.1,.2,.3,.4)
#R <- diag(4)
#Rsq <- t(Beta) %*% R %*% Beta
#Rsq
#out<-fungibleR(R,Beta,Lp = .12345, eps=1e-16)
#round(out$RstarLp,5)
#eigen(out$RstarLp)$values



# R <- matrix(.5,5,5)
# diag(R) <- 1
# Beta <- rep(.1,5)
# # 
# # #change value of Lp to control the size of the smallest eigenvalue of RstarMax
# out<-fungibleR(R,Beta,Lp = 0.12345678, eps=1e-16)
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


