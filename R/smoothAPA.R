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



