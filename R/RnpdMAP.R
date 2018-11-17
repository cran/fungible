RnpdMAP <- function (Rpop, 
                     Lp = NULL, 
                     NNegEigs = 1,
                     NSmoothPosEigs = 4,
                     NSubjects = NULL, 
                     NSamples = 0, 
                     MaxIts = 15000, 
                     PRINT=FALSE, 
                     Seed = NULL) 
{
  #############################################################################
  ## Sample NPD R matrices from a user-supplied
  ## population R using the alternating projections algorithm
  ##
  ## Niels Waller
  ## Dec 12 2017 3:10pm
  ##
  ## Arguments:
  ##  Rpop              Input (PD or PSD) p x p Population correlation matrix
  ##
  ##  Lp                Desired minimum eigenvalue in the NPD matrices
  ##
  ##  NNegEigs          Number of eigenvalues < 0 in Rnpd
  ##
  ##
  ##  NSmoothPosEigs    Number of eigenvalues > 0 to smooth:  the smallest 
  ##                    NSmoothPosEigs > 0  be smoothed toward 0.
  ##
  ##  NSubjects         Sample size (required when NSamples > 0) parameter used to 
  ##                    generate sample correlation matrices. Default = NULL.
  ##
  ##  NSamples          Generate NSamples sample R matrices. If NSamples = 0
  ##                    the program will attempt to find Rnpd such that 
  ##                    ||Rpop - Rnpd||_2 is minimized. 
  ##
  ##  MaxIts            maximum number of projection interations.
  ##
  ##  PRINT             (Logical) If TURE the program will print the iteration history for Lp. 
  ##                    Default = NULL. 
  ##  Seed              Optional seed for random number generation.
  ##
  ##  Values:
  ##  Rpop    Population correlation matrix.
  ##
  ##  R       Sample correlation matrix.
  ##
  ##  Rnpd    NPD improper (pseudo) correlation matrix.
  ##
  ##  Lp      desired value of minimum eigenvalue.
  ##
  ##  minEig  observed value of minimum eigenvalue of Rnpd.
  ##
  ##  convergence     0 = converged
  ##                  1 = not converged in MaxIts.
  ##
  ##  feasible       (Logical) TRUE if max(abs(r_ij)) <=1  
  ##
  ##  Seed          saved seed for random number generator.
  ##
  ##  prbs1
  ##
  ##  prbs2
  ##########################################################################
  
  
  if(NSmoothPosEigs == 0)  NSmoothPosEigs <- NULL
  prbs2 <- NULL
  
  if(is.null(Seed)){
    Seed<- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) 
  }
  set.seed(Seed)
  
  ## Generate monotonically decreasing negative eigenvalues with 
  ## min eig = Lp
  prbs1Fnc <- function(){
    prbs <- 1
    if(NNegEigs > 1){
          prbs <- 1- c(log(sort(runif(NNegEigs-1,1,100), decreasing = TRUE))/log(100),0)
    } #end if else
    prbs
   } ## END prbsFnc
  
  
  if( NSamples > 0 ){
   ## If NSamples > 0 then Generate sample correlation matrices from Wishart distribution
   WishartArray <- rWishart(n = NSamples, df = NSubjects-1, Sigma = Rpop)
   WishartList <-lapply( seq(dim(WishartArray)[3]), function(x) WishartArray[ , , x] )
   ## scale to corr matrices with cov2cor
   corList<-lapply(WishartList,cov2cor)
  }
  else if (NSamples == 0)
   corList <- list(Rpop)
  
  
## Lpop = eigenvalues of population PD R
  Lpop <- eigen(Rpop, only.values = TRUE)$values
  
 
 
## Project matrix onto Symmetric A with min eig Lp 
## and (NNegEigs) smoothed eigs < 0  
  proj.S.Lp <- function(A, Lp, Nvar, prbs1, prbs2){

    VLV <- eigen(A)
    V <- VLV$vectors
    L <- VLV$values
    L[L < Lp] <- Lp
    
    if(NNegEigs == 1){
     L[Nvar] <- Lp
     } else{
      #  Smooth the smaller eigevalues < 0
      L[ (Nvar - (NNegEigs - 1)):Nvar] <- prbs1 * Lp
      
      # Smooth the smaller eigevalues > 0
      if(!is.null(NSmoothPosEigs)){
        L[(Nvar - (NNegEigs + NSmoothPosEigs )):(Nvar - (NNegEigs + 1 ))] <- 
           prbs2 * Lpop[(Nvar - (NNegEigs + NSmoothPosEigs )):(Nvar - (NNegEigs + 1 ))]
      }
    } # End else
   
    # reconstruct symmetric Cov matrix
    A <- V %*% diag(L)%*% t(V)
    (A + t(A))/2
  } #End  proj.S.Lp
  
## Project onto symetric matrix set with unit diagonals in bounding cube 
  proj.U <- function(A){
    diag(A) <- 1
    
    # find row col indices for |rij| > 1
    # and bring elements into cube
    rc.ind <- which(abs(A) > 1, arr.ind = TRUE)
    A[rc.ind] <- sign(A[rc.ind]) 
    A
    
  }

  
## Find NPD pseudo R matrix
   RAPA <- function(R, Lp, MaxIts = 1000){
    
      Nvar <- ncol(R)
  
      EigTest <- FALSE
      tries <- 0
      converged <- TRUE 
      EigHx <- rep(0, MaxIts)
  
      Rk <- R
    
     # prbs1 and prbs2 are use to smooth the rate of 
     # descrease in the smallest eigenvalues.
      
     # prbs1 smooths NNegEigs eigs < 0
     prbs1 <- prbs1Fnc()
    
     # prbs2 smooths the smallest NSmoothPosEigs eigs > 0
     if(!is.null(NSmoothPosEigs)){
       prbs2 <- log(sort(runif(NSmoothPosEigs,1,100), decreasing = TRUE))/log(100)
     }
     
     while (!EigTest) 
     {
        tries <- tries + 1
        if (tries > MaxIts) {
           converged = FALSE
           break
        }
      
     ## Project onto S with min Eig Lp
     Sk <- proj.S.Lp(Rk, Lp, Nvar, prbs1, prbs2) 
    
     ## project onto space of Symmetric matrices, Uk, with unit diagonal (in bounding hypercube) 
     Uk <- proj.U(Sk)
     
     # Rk+1 = Uk
     Rk <- Uk
    
     minEig <- eigen(Rk, symmetric = TRUE, only.values = TRUE)$val[Nvar]
     EigHx[tries] <- minEig
     
     if(PRINT){
       cat("\nAt iter ", tries, " min eig = ", eigen(Rk)$val[Nvar])
     }
     
     ## test convergence
     if(1e+10*(minEig - Lp)^2 < 1e-12)   #
       EigTest <- TRUE
     } ## END while loop
  
     EigHx <- EigHx[1:tries]
     colnames(Rk) <- rownames(Rk) <- colnames(R)
  
     feasible <- TRUE
  ## if any |rxy| > 1 then reject solution 
     if(max(abs(Rk)) > 1) feasible <- FALSE

     list(Rpop = Rpop, 
         R=R, 
         Rnpd = Rk, 
         Lp = Lp, 
         ObsLp = minEig, 
         EigHx = EigHx, 
         converged = converged,
         feasible = feasible,
         Seed = Seed,
         prbs1 = prbs1,
         prbs2 = prbs2)
     } #End internal function RAPA
  
  ## Call Rnpd and return output as a list
  lapply(corList, RAPA, Lp, MaxIts)
}  ## END FUNCTION


#out <- RnpdMAP(Rpop=Rpop, Lp = -.5, NSubjects = 200, NSamples = 10, MaxIts = 5000, PRINT=FALSE) 




