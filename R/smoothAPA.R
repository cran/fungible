#' Smooth a NPD R matrix to PD using the Alternating Projection Algorithm
#' 
#' Smooth a Non positive defnite (NPD) correlation matrix to PD using the
#' Alternating Projection Algorithm with Dykstra's correction via Theory
#' described in Higham 2002.
#' 
#' 
#' @param R A p x p indefinite matrix.
#' @param delta Desired value of the smallest eigenvalue of smoothed matrix,
#' RAPA. (Default = 1e-06).
#' @param fixR User-supplied integer list that instructs the program to
#' constrain elements in RAPA to equal corresponding elements in RAPA. For
#' example if fixR = c(1,2) then smoothed matrix, RAPA[1:2,1:2] = R[1:2,1:2].
#' Default (fixR = NULL).
#' @param Wghts A p-length vector of weights for differential variable
#' weighting. Default (Wghts = NULL).
#' @param maxTries Maximum number of iterations in the alternating projections
#' algorithm. Default (maxTries = 1000).
#' @return \item{RAPA}{A smoothed matrix.} \item{delta}{User-supplied delta
#' value.} \item{Wghts}{User-supplied weight vector.} \item{fixR}{User-supplied
#' integer list that instructs the program to constrain elements in RAPA to
#' equal corresponding elements in R.} \item{convergence}{A value of 0
#' indicates that the algorithm located a feasible solution. A value of 1
#' indicates that no feasible solution was located within maxTries.}
#' @author Niels Waller
#' @keywords statistics
#' @export
#' @examples
#' 
#' data(BadRKtB)
#' 
#' ###################################################################
#' ##  Replicate analyses in Table 2 of Knol and ten Berge (1989).
#' ###################################################################
#' 
#' ## n1 = 0,1
#' out<-smoothAPA(R = BadRKtB, delta = .0, fixR = NULL, Wghts = NULL, maxTries=1e06)
#' S <- out$RAPA
#' round(S - BadRKtB,3)
#' normF(S - BadRKtB)
#' eigen(S)$val
#' 
#' ## n1 = 2
#' out<-smoothAPA(R = BadRKtB, fixR =c(1,2), delta=.0, Wghts = NULL, maxTries=1e06)
#' S <- out$RAPA
#' round(S - BadRKtB,3)
#' normF(S - BadRKtB)
#' eigen(S)$val
#' 
#' ## n1 = 4
#' out<-smoothAPA(R = BadRKtB, fixR = 1:4, delta=.0, Wghts = NULL, maxTries=1e06)
#' S <- out$RAPA
#' round(S - BadRKtB,3)
#' normF(S - BadRKtB)
#' eigen(S)$val
#' 
#' ## n1 = 5
#' out<-smoothAPA(R = BadRKtB, fixR = 1:5, delta=0, Wghts = NULL, maxTries=1e06)
#' S <- out$RAPA
#' round(S - BadRKtB,3)
#' normF(S - BadRKtB)
#' eigen(S)$val
#' 
#' ###################################################################
#' ##  Replicate analyses in Table 3 of Knol and ten Berge (1989).
#' ###################################################################
#' 
#' ## n1 = 0,1
#' out<-smoothAPA(R = BadRKtB, delta = .05, fixR = NULL, Wghts = NULL, maxTries=1e06)
#' S <- out$RAPA
#' round(S - BadRKtB,3)
#' normF(S - BadRKtB)
#' eigen(S)$val
#' 
#' ## n1 = 2
#' out<-smoothAPA(R = BadRKtB, fixR =c(1,2), delta=.05, Wghts = NULL, maxTries=1e06)
#' S <- out$RAPA
#' round(S - BadRKtB,3)
#' normF(S - BadRKtB)
#' eigen(S)$val
#' 
#' ## n1 = 4
#' out<-smoothAPA(R = BadRKtB, fixR = 1:4, delta=.05, Wghts = NULL, maxTries=1e06)
#' S <- out$RAPA
#' round(S - BadRKtB,3)
#' normF(S - BadRKtB)
#' eigen(S)$val
#' 
#' ## n1 = 5
#' out<-smoothAPA(R = BadRKtB, fixR = 1:5, delta=.05, Wghts = NULL, maxTries=1e06)
#' S <- out$RAPA
#' round(S - BadRKtB,3)
#' normF(S - BadRKtB)
#' eigen(S)$val
#' 
#' ###################################################################
#' ## This example illustrates differential variable weighting.
#' ## 
#' ## Imagine a scenerio in which variables 1 & 2 were collected with 
#' ## 5 times more subjects than variables 4 - 6 then . . .
#' ###################################################################
#' ## n1 = 2
#' out<-smoothAPA(R = BadRKtB, delta=.0, fixR = NULL, Wghts = c(5, 5, rep(1,4)), maxTries=1e5)
#' S <- out$RAPA
#' round(S - BadRKtB,3)
#' normF(S - BadRKtB)
#' eigen(S)$val
#' 
smoothAPA <- function(R, delta = 1e-06, fixR = NULL, Wghts = NULL, maxTries = 1000){
  ## Modified Alternating Projection algorithm with Dykstra's 
  ## correction based on theory described in:
  ## Higham, N. J. (2002). Computing the nearest correlation matrix: 
  ## A problem from finance. IMA Journal of Numerical Analysis, 
  ## 22(3), 329--343. 
  ##
  ##  Arguments
  ##   R         :   A p x p improper correlation matrix.
  ##   delta     :   Desired value of the smallest eigenvalue of smoothed matrix, RAPA. (Default = 1e-06). 
  ##   fixR      :   A vector of integers that specifies which variables to hold fixed. For example
  ##                if fixR = c(1,2)  then smoothed matrix RAPA[1:2,1:2] = R[1:2,1:2]. 
  ##                Default (fixR = NULL).
  ##   maxTries  :  Maximum number of iterations in the alternating projections algorithm.
  ##   Wghts     :   A p-length vector of weights for differential variable weighting. Default (Wghts = NULL).
  ##
  ##  Value
  ##   RAPA      :   A smoothed correlation matrix
  ##   delta     :   User-supplied value of delta
  ##   Wghts     :   User-supplied vector of weights
  ##   fixR      :   User-supplied vector of integers specifying fixed correlations 
  
  
  if(min(eigen(R)$val) > 0 ) stop("\nFATAL ERROR:  No smoothing necessary: Input matrix is not indefinite\n")
    if(!is.null(Wghts) & abs(delta) > 0){
     stop("\nFATAL ERROR: delta must equal 0 when applying weights (Wghts)\n", call.=FALSE)
    }
  
  Ck <- R
  Nvar <- ncol(R)
  DeltaS <- matrix(0,nrow=Nvar,ncol=Nvar)
  
  
  ## Differentially weight rij
  if(!is.null(Wghts) & length(Wghts) != Nvar) stop("\nWeight vector of wrong length\n")
  if(!is.null(Wghts)) {
     ## check that all weights are positive
     if(min(Wghts) <= 0) stop("\nAll weights must be positive\n")
     Dw <- diag(sqrt(Wghts))
     DwInv <- diag(1/sqrt(Wghts))
    }  
  else
    Dw <-  DwInv <-  diag(Nvar)
   

  
  NPD <- TRUE
  convergence <- -99
  tries <- 0
###########################  
#########  BEGIN WHILE LOOP   
  while(NPD & (tries <= maxTries)){
    
    tries <- tries + 1
    if(tries >= maxTries) {
      print(eigen(Rk)$val)
      convergence <- 1
      warning("exceeded maximum tries")

    }
    #Dykstra's correction
    Rk <- Ck - DeltaS
  
    
    
    # compute eigen decomposition of possibly weighted matrix
    VLV <- eigen( Dw %*% Rk %*% Dw)
    V <- VLV$vec
    L<- VLV$val
    L[L<=0] <- delta
    Lplus <- diag(L)
    
    
    # Project into space of PSD or PD symmetric matrices
    Ck <- DwInv %*% V %*% Lplus %*% t(V) %*% DwInv
    
    DeltaS <- Ck - Rk
    
    # Project into space of symmetric matrices with unit diagonal values
    diag(Ck) <- 1
	

    ## Fix elements of RAPA to corresponding elements of R
    if(!is.null(fixR)) Ck[fixR,fixR] <- R[fixR,fixR]
    
  
		# If smallest eigenvalue == delta (within tolerance) End while loop
     if( 5e05*(eigen(Ck)$val[Nvar] - delta )^2 < 1e-12) {
       NPD <- FALSE
       convergence <- 0 # problem converged
     } 
    
  }
  
  ## Reconstruct smoothed matrix at solution
  Ck <- DwInv %*% V %*% Lplus %*% t(V) %*% DwInv
  colnames(Ck) <- rownames(Ck) <- colnames(R)
  
  if(tries >= maxTries)  convergence <- 1
  
  list(RAPA = Ck, delta = delta, Wghts=Wghts, fixR = fixR, convergence = convergence)
}  



