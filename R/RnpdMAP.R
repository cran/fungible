#' Generate Random NPD R matrices from a user-supplied population R
#' 
#' Generate a list of Random NPD (pseudo) R matrices with a user-defined fixed
#' minimum eigenvalue from a user-supplied population R using the method of
#' alternating projections.
#' 
#' 
#' @param Rpop input (PD or PSD) p x p Population correlation matrix.
#' @param Lp desired minimum eigenvalue in the NPD matrices.
#' @param NNegEigs number of eigenvalues < 0 in Rnpd.
#' @param NSmoothPosEigs number of eigenvalues > 0 to smooth: the smallest
#' NSmoothPosEigs > 0be smoothed toward 0.
#' @param NSubjects sample size (required when NSamples > 0) parameter used to
#' generate sample correlation matrices. Default = NULL.
#' @param NSamples generate NSamples sample R matrices. If NSamples = 0 the
#' program will attempt to find Rnpd such that ||Rpop - Rnpd||_2 is minimized.
#' @param MaxIts maximum number of projection iterations.
#' @param PRINT (logical) If TURE the program will print the iteration history
#' for Lp. Default = NULL.
#' @param Seed Optional seed for random number generation.
#' @return \item{Rpop}{population (PD) correlation matrix.} \item{R}{sample
#' correlation matrix.} \item{Rnpd}{NPD improper (pseudo) correlation matrix.}
#' \item{Lp}{desired value of minimum eigenvalue.} \item{minEig}{observed value
#' of minimum eigenvalue of Rnpd.} \item{convergence}{0 = converged; 1 = not
#' converged in MaxIts iterations of the alternating projections algorithm.}
#' \item{feasible}{logical) TRUE if max(abs(r_ij)) <= 1. If FALSE then one or
#' more values in Rnpd > 1 in absolute value.} \item{Seed}{saved seed for
#' random number generator.} \item{prbs1}{vector probabilities used to generate
#' eigenvalues < 0.} \item{prbs2}{vector of probabilities used to smooth the
#' smallest NSmoothPosEigs towards zero.}
#' @author Niels G. Waller
#' @keywords fungible
#' @export
#' @import MASS
#' @examples
#' 
#' library(MASS)
#' 
#' Nvar = 20
#' Nfac = 4
#' NSubj = 600
#' Seed = 123    
#' 
#' set.seed(Seed)
#' 
#' ## Generate a vector of classical item difficulties
#' p <- runif(Nvar)
#' 
#' cat("\nClassical Item Difficulties:\n")
#' 
#' print(rbind(1:Nvar,round(p,2)) )
#' 
#' summary(p)
#' 
#' 
#' ## Convert item difficulties to quantiles
#' b <- qnorm(p)
#' 
#' 
#' ## fnc to compute root mean squared standard deviation
#' RMSD <- function(A, B){
#'   sqrt(mean( ( A[lower.tri(A, diag = FALSE)] - B[lower.tri(B, diag = FALSE)] )^2))
#' }
#' 
#' 
#' ## Generate vector of eigenvalues with clear factor structure
#'   L <- eigGen(nDimensions = Nvar, 
#'             nMajorFactors = Nfac, 
#'             PrcntMajor = .60, 
#'             threshold  = .50)
#'           
#' 
#' ## Generate a population R matrix with the eigenvalues in L
#'   Rpop <- rGivens(eigs = L)$R
#'   
#' ## Generate continuous data that will reproduce Rpop (exactly)
#'   X <- mvrnorm(n = NSubj, mu = rep(0, Nvar), 
#'                Sigma = Rpop, empirical = TRUE)
#'                
#' if( any(colSums(X) == 0) ){
#'   stop("One or more variables have zero variance. Generate a new data set.")               
#'  }
#'  
#' ## Cut X at thresholds given in b to produce binary data U
#'   U <- matrix(0, nrow(X), ncol(X))
#'   for(j in 1:Nvar){
#'     U[X[,j] <= b[j],j] <- 1
#'   }
#'   
#' ## Compute tetrachoric correlations
#'   Rtet <- tetcor(U, Smooth = FALSE, PRINT = TRUE)$r
#'   # Calculate eigenvalues of tetrachoric R matrix
#'   Ltet <- eigen(Rtet)$values
#'   
#'   if(Ltet[Nvar] >= 0) stop("Rtet is P(S)D")
#'   
#' ## Simulate NPD R matrix with minimum eigenvalue equal to 
#'   # min(Ltet)
#'   out <- RnpdMAP(Rpop, 
#'                Lp = Ltet[Nvar], 
#'                NNegEigs = Nvar/5,
#'                NSmoothPosEigs = Nvar/5, 
#'                NSubjects = 150, 
#'                NSamples = 1, 
#'                MaxIts = 15000, 
#'                PRINT = FALSE, 
#'                Seed = Seed) 
#' 
#' ## RLp is a NPD pseudo R matrix with min eigenvalue = min(Ltet)
#'   RLp <- out[[1]]$Rnpd
#' 
#' ## Calculate eigenvalues of simulated NPD R matrix (Rnpd)
#'   Lnpd <- eigen(RLp, only.values = TRUE)$values
#'   
#' ## Scree plots for observed and simulated NPD R matrices.  
#'   ytop <- max(c(L,Lnpd,Ltet))
#'   pointSize = .8
#'   plot(1:Nvar, L, typ = "b", col = "darkgrey", lwd=3, 
#'        lty=1, 
#'        main = 
#'        "Eigenvalues of Rpop, Tet R, and Sim Tet R:
#'        \nSimulated vs Observed npd Tetrachoric R Matrices",
#'        ylim = c(-1, ytop),
#'        xlab = "Dimensions", 
#'        ylab = "Eigenvalues",
#'        cex = pointSize,cex.main = 1.2)
#'   points(1:Nvar, Lnpd, typ="b", 
#'          col = "red", lwd = 3, lty=2, cex=pointSize)
#'   points(1:Nvar, Ltet, typ="b", 
#'          col = "darkgreen", lwd = 3, lty = 3, cex= pointSize)
#'  
#'   legend("topright", 
#'          legend = c("eigs Rpop", "eigs Sim Rnpd", "eigs Emp Rnpd"), 
#'          col = c("darkgrey", "red","darkgreen"), 
#'          lty = c(1,2,3), 
#'          lwd = c(4,4,4), cex = 1.5)
#'   
#'   abline(h = 0, col = "grey", lty = 2, lwd = 4)
#'  
#'   cat("\nRMSD(Rpop, Rtet) = ", round(rmsd(Rpop, Rtet), 3))
#'   cat("\nRMSD(Rpop, RLp) = ",  round(rmsd(Rpop, RLp),  3))
#' 
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




