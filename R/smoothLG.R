

smoothLG <- function(R, start.val = NULL, Wghts = NULL, PD = FALSE, 
	        Penalty = 50000, eps=1e-07){
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
## smoothLG: a function for smoothing an indefinite  matrix 
## to a PSD matrix via theory described by Lurie and Goldberg
## Niels Waller
## February 14, 2016
## February 29, 2016 added var names
## June 30, 2016 added error checks
##
## Arguments:
## R      Indefinite Matrix
## start.val  optional start values for Cholesky factor of S  
## Wghts     an optional matrix of weights such that
##           fnc minimizes wij(rij - sij)^2
##           where wij is Wght[i,j]  
## PD        logical. if PD = TRUE the fnc will smooth
##           the least squares solution to insure PD
## Penality  Scalar weight for the Lagrangian multiplier.
##           Defaut = 50000
## eps       value to add to zero eigenvalues if smoothed matrix must be PD  

##
## Value:
##        RLG           Lurie Goldberg smooth matrix
##        RKB           Knol and Berger smoothed matrix
##        convergence   0 = converged, 1 = failure
##        start.val     start.values
##        gr            analytic gradient at solution 
##        Penalty  
##        PD            user-supplied value of PD  
##        Wghts
##        eps  
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
  
## add some error checkig code to check that input 
##  matrices are symmetric etc.  
# if(!isSymmetric(R)) stop("\nInput matrix not symmetric\n") 
  
if( (PD!=TRUE) & (PD!=FALSE) ) stop("\nPD should be a logical (TRUE or FALSE)\n") 
  
## Declarations and definitions
LG = 1  #start val for Lagrange multiplier    
Nvar <- ncol(R)  
lowtriN<-Nvar*(Nvar+1)/2
I <- diag(Nvar)


## if weights not supplied all weights are 1
if(is.null(Wghts)) W <-matrix(1,Nvar,Nvar)
else W <- Wghts

## Z is a matrix of Zeros
Z <-matrix(0,Nvar,Nvar)
## define trace function 
TR <- function(x) sum(diag(x))
## define custom diag function
DIAG <- function(x) diag(diag(x))  


##  corSmooth: a function to smooth an ID matrix to a PD matrix 
## by a function described by Knol and Berger
## This function is used to generate start values
corSmooth <-function (R, eps = .01)    #1e+08 * .Machine$double.eps) 
{
  KDK <- eigen(R)
  Dplus <- D <- KDK$values
  Dplus[D <= 0] <- eps
  Dplus <- diag(Dplus)
  K <- KDK$vectors
  Rtmp <- K %*% Dplus %*% t(K)
  invSqrtDiagRtmp <- diag(1/sqrt(diag(Rtmp)))
  invSqrtDiagRtmp %*% Rtmp %*% invSqrtDiagRtmp
}

## smooth R by Knol and Berger method
## Rkb = R Knol and Berger
RKB<-corSmooth(R)

## Compute Cholesky of RS for start values
if( is.null(start.val) ){
   start.val <-c((t(chol(RKB)))[lower.tri(RKB, diag=T)][-1], LG)
}


## position of the Lagrange multiplier in the
## list of parameter estimates
LagrangeNum <- length(start.val)
 

## makeL create a lower-triangular matrix of the same 
## dimension as R. This will be updated during optimization.
  makeL <- function(x){
    L <- I
    lower <- lower.tri(L, diag = TRUE)
    L[lower] <- c(1,x[-LagrangeNum])
    L
   } #End makeL

## fnc to minimize (1/2) squared euclidean distance
## between S and  R. 
  
  fnc <- function(x){
    L <- makeL(x)
    S <- L%*% t(L)
    .5*TR( (sqrt(W) * (R - S)) %*% ( sqrt(W) *(R - S)) ) + x[LagrangeNum]* Penalty* TR( DIAG(S-I) %*% DIAG(S-I) ) 
  } #End fnc
  

###############################################
## Compute gradient of fnc
grvec<-rep(0,lowtriN-1)

grFNC <- function(x){
  
  L <- makeL(x)
  Lt <- t(L)
  S <- L %*% t(L)
  
## function to compute gradient element Lij
  gradLij <- function(i,j){
  
    Jij<-Jji<-Z
    Jij[i,j]<-1
    Jji[j,i]<-1

## partial deriv of main function w.r.t. Lij   
    dfdLij <- TR( (W*(R - S)) %*% -(L%*%Jji + Jij%*%Lt))
## partial deriv of constraint function w.r.t. Lij   
    dgdLij <- 2 * x[lowtriN]* Penalty * TR( (DIAG(S - I) %*% (L%*%Jji + Jij%*%Lt)) )   

    dfdLij + dgdLij
} #End gradLij

## fill in gradient vector
   k<-0
    for(j in 1:Nvar){
      for(i in 2:Nvar){
        if(i>=j) {
          k = k+1
          grvec[k]<-gradLij(i,j)
         }   
      }
    } 
   grvec[k+1] <- Penalty * TR( DIAG(S-I)%*%DIAG(S-I) )
   grvec
} 


##minimize fnc.
  out <- optim(par = start.val, fn = fnc, 
               gr=grFNC,
               method="BFGS", 
               control=list(fnscale=1,
                            maxit=10000,
                            reltol=1e-20, 
                            abstol=1e-20))
 
  ## check if function converged 

  if(normF(grFNC(out$par)) > 1e-05) out$convergence <-1
 
  
  # assembled smoothed matrix
    L<-makeL(out$par[-LagrangeNum])
    Rsmooth <- L %*% t(L)
    
    
    
    ## Is smoothed matrix a proper correlation matrix
    if( abs(sum(abs(diag(Rsmooth))) - Nvar) > 1e-05 ) out$convergence <-1
    if(det(Rsmooth) > 1) out$convergence <-1
  
  if(out$convergence==0){  
    
    ## Smooth to insure PD matrix?
    if(PD) Rsmooth<-corSmooth(Rsmooth, eps)


    colnames(Rsmooth) <- rownames(Rsmooth) <- colnames(R)
    colnames(RKB) <- rownames(RKB) <- colnames(R)
  
    list(RLG = Rsmooth, 
         RKB = RKB, 
         convergence = out$convergence,
         par = out$par,
         start.val = start.val,
         gr=grFNC(out$par),
         Penalty = Penalty,
         PD = PD,
         Wghts = W,
         eps = eps
         )
  }
  ## if function failed to converge 
   else list(RLG = 999, 
             RKB = 999, 
             convergence = 1,
            # par = rep(99, length(start.val)),
             start.val = start.val,
             gr=NULL,
             Penalty = Penalty,
             PD = PD,
             Wghts=W,
             eps=eps
             )
 
} ## End smoothLG


