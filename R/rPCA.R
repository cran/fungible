## Niels Waller
## rPCA
## June 6, 2021
## Sept 21, 2020
## Generate a population R matrix for a 
## user-defined PCA model via the 
## alternating projections algorithm

#' Generate a Correlation Matrix from a Truncated PCA Loadings Matrix. 
#' 
#' This function generates a random (or possibly unique) correlation matrix (R) from an unrotated 
#' or orthogonally rotated PCA loadings matrix via a modified alternating 
#' projections algorithm.
#' 
#' @param F (Matrix) A p (variables) by k (components) PCA loadings matrix. 
#'    F can equal either an unrotated or an orthogonally rotated 
#'    loadings matrix.
#' @param epsMax (Scalar) A small number used to evaluate function convergence. 
#'     Default (epsMax = 1E-18).
#' @param maxit (Integer) An integer that specifies the maximum  number of 
#' iterations of the modified alternating projections algorithm (APA).
#' @param  Seed (Integer) A user-defined starting seed for the 
#' random number generator.  If Seed = NULL then rPCA will generate a 
#' random starting seed. Setting Seed to a positive integer will generate 
#' reproducible results.  Default (Seed = NULL) 
#' @param InitP2 (Integer) The method used to initiate the remaining columns of the
#' truncated principal components solution.  If \code{InitP2 = 1} then the 
#' starting P2 will be a random semi-orthogonal matrix.  If If \code{InitP2 = 2} then 
#' the starting P2 will be a semi-orthogonal matrix that is in the left 
#' null space of P1. Default (\code{InitP2 = 2}). Of the two options, 
#' \code{InitP2 = 2} generally converges to a single feasible solution in less time.
#' \code{InitP2 = 1} can be used to generate different solutions from different 
#' starting seeds. 
#' @param Eigs (Vector) Under some conditions, \code{rPCA} can generate (or 
#' reproduce)  a unique correlation matrix with known (i.e., user-specified) 
#' eigenvalues from a truncated PC loadings matrix, \code{F}, even when 
#' the rank of \code{F} is less than p (the number of observed variables). 
#' \code{Eigs} is an optional p-length vector of eigenvalues for \code{R}. 
#' Default (\code{Eigs = NULL}).
#' @param PrintLevel (Integer) If PrintLevel = 0 no output will be printed
#'   (choose this option for Monte Carlo simulations). If PrintLevel = 1 the 
#'   program will print the APA convergence status and the number of iterations 
#'   used to achieve convergence.  If PrintLevel = 2 then rPCA will print the 
#'   iteration convergence history of the modified APA algorithm. Default 
#'   (PrintLevel = 1).
#'   
#' @return 
#' \itemize{
#'   \item \strong{R} (Matrix) A p by p correlation matrix that 
#'    generates the desired PCA loadings.
#'  \item \strong{Tmat} (Matrix) A k by k orthogonal rotation matrix 
#'    that will rotate the unrotated PCA loadings matrix, P1, to F (if F is an 
#'    orthogonally rotated loadings matrix).
#'   \item \strong{P1} (Matrix)  The p by k unrotated PCA loadings matrix 
#'    that is associated  with F.
#'  \item \strong{Fhat} (Matrix) The p by k estimated (and possibly rotated) 
#'   PCA loadings matrix from the simulated matrix R.
#'  \item \strong{error}  (Logical) A logical that indicates whether 
#'    F is a legitimate PCA loadings matrix. 
#'  \item \strong{Lambda} (Vector) The sorted eigenvalues of R.
#' \item \strong{iterHx} (Vector) Criterion (i.e., fit) values for  
#'  for each iteration of the modified APA algorithm.
#' \item \strong{converged} (Logical) A logical that signifies function 
#'    convergence.
#' \item \strong{Seed} (Integer) Either a user-defined or function generated 
#'   starting seed for the random number generator. 
#'  }
#'  
#' @references Escalante, R. and Raydan, M.  (2011).  Alternating projection 
#'  methods.   Society for Industrial and Applied Mathematics. 
#' @references ten Berge, J. M. and Kiers, H. A.  (1999).  Retrieving the 
#' correlation matrix from a truncated PCA solution: The inverse principal 
#' component problem.  Psychometrika, 64(3), 317--324.  

#' @author Niels G. Waller (nwaller@umn.edu)
#' @examples
#' 
#' # External PCA function ---
#' # used to check results
#'  
#' PCA <- function(R, k = NULL){
#'   if(is.null(k)) k <- ncol(R)
#'   VLV <- eigen(R)
#'  V <- VLV$vectors
#'  L <- VLV$values
#' 
#'  
#'  if( k > 1){
#'    P <-  V[, 1:k] %*% diag(L[1:k]^.5)
#'  }
#'  else{
#'    P <- as.matrix(V[, 1], drop=False) * L[1]^.5
#'  }
#'   Psign <- sign(apply(P, 2, sum))
#'   if(k > 1) Psign = diag(Psign)
#'   P <- P %*%  Psign
#'  P
#' }#END PCA  
#' 
#'   
#' ## Generate Desired Population rotated PCA loadings matrix
#' ## Example = 1
#'  k = 2
#'  F <- matrix(0, 8, 2) 
#'  F[1:4, 1] <- seq(.75, .72, length= 4)  
#'  F[5:8, 2] <- seq(.65, .62, length= 4)  
#'  F[1,2] <- .1234
#'  F[8,1] <- .4321
#'  colnames(F) <-   paste0("F", 1:k) 
#'  (F)
#'  
#'  ## Run Example 1
#'  pout <- rPCA(F, 
#'               maxit = 5000, 
#'               Seed = 1, 
#'               epsMax = 1E-18,
#'               PrintLevel = 1)
#' pout$converged
#' eigen(pout$R)$values
#' if(pout$error == FALSE & pout$converged){ 
#'     Fhat <- pout$Fhat
#'     cat("\nPCA Loadings\n")
#'     ( round( cbind(F,Fhat ), 5) )
#'  }
#'  
#'  ## Example = 2      
#'  ## Single component example from Widaman 2018
#'
#'  k = 1
#'  F <- matrix(rep(c(.8,.6, .4), each = 3 ), nrow = 9, ncol = 1)
#'  colnames(F) <-   paste0("F", 1:k) 
#'  (F)
#' 
#'  ## Run Example 2
#'  pout <- rPCA(F, 
#'               maxit = 5000, 
#'               Seed = 1, 
#'               epsMax = 1E-18,
#'               PrintLevel = 1)
#'  pout$converged
#'  pout$Fhat
#'  eigen(pout$R)$values
#' if(pout$error == FALSE & pout$converged){ 
#'     Fhat <- pout$Fhat
#'     cat("\nPCA Loadings\n")
#'     ( round( cbind(F,Fhat ), 5) )
#'  }
#'   
#' ## Example 3 ----
#' ## 2 Component example from Goldberg and Velicer (2006).
#'  k = 2
#'  F = matrix(c( .18, .75,
#'                .65, .19,
#'                .12, .69,
#'                .74, .06,
#'                .19, .80,
#'                .80, .14,
#'               -.05, .65,
#'                .71, .02), 8, 2, byrow=TRUE)
#'  colnames(F) <-   paste0("F", 1:k) 
#'  (F)
#'
#' ## Run Example 3
#' pout <- rPCA(F, 
#'             maxit = 5000, 
#'             Seed = 1, 
#'             epsMax = 1E-18,
#'             PrintLevel = 1)
#' pout$converged
#'
#' eigen(pout$R)$values
#'
#' if(pout$error == FALSE & pout$converged){ 
#'   Fhat <- pout$Fhat
#'   cat("\nPCA Loadings\n")
#'   ( round( cbind(F,Fhat ), 5) )
#' #
#' #
#' ## Example 4
#' #
#' SEED = 4321
#' set.seed(SEED)
#' k= 3
#' ## Generate eigenvalues for example R matrix
#' L7 <- eigGen(nDimensions = 7,
#'              nMaj = 3,
#'              PrcntMajor = .85,
#'              threshold = .8)
#' 
#' ## Scree Plot
#' plot(1:7, L7, 
#'     type = "b", 
#'     ylim = c(0,4),
#'     main = "Scree Plot for R",
#'     ylab = "Eigenvalues",
#'     xlab = "Dimensions")
#'
#' ## Generate R
#' R <- rGivens(eigs=L7, Seed = SEED)$R
#' print( R, digits = 4)
#'
#' #Extract loadings for 3 principal components
#' F <- PCA(R, k = k)
#' 
#' # rotate loadings with varimax to examine underlying structure
#' print( round(varimax(F)$loadings[], 3) )
#' 
#' ## run rPCA with user-defined eigenvalues
#' rout <- rPCA(F,
#'             epsMax = 1e-20, 
#'             maxit = 25000, 
#'             Seed = SEED,   
#'             InitP2 = 1,
#'             Eigs = L7,
#'             PrintLevel = 1) 
#'
#' ## Compute PCA on generated R
#' 
#' Fhat <- PCA(rout$R, k = 3)
#' #
#' ## align factors
#' Fhat <- fungible::faAlign(F, Fhat)$F2
#' 
#' ## Compare solutions
#' print( round( cbind(F, Fhat), 5) )
#'
#' ## Compare Eigenvalues
#' print( cbind(L7, eigen(rout$R)$values ), digits=8) 
#' #
#' ## Compare R matrices: 8 digit accuracy
#' print( round(R - rout$R, 8) )
#'
#' }
#' @export

# rPCA ----
rPCA <- function(F, 
                 epsMax = 1E-18,  
                 maxit = 2000,   
                 Seed = NULL, 
                 InitP2 = 2,     
                 Eigs = NULL,
                 PrintLevel = 1){
  
  
    # Dimensions of F
    p <- nrow(F)
    k <- ncol(F)
    
    if(is.null(Eigs)) EigsSqrt = NULL
    # sqrt of the eigs in the last (p-k) components
    if(!is.null(Eigs)) EigsSqrt = sqrt(Eigs[(k+1):p])
    
    # Generate random seed if not supplied
    if(is.null(Seed)) Seed <- sample(1:10000, 1)
  
    # iterHx will contain iteration history
    iterHx <- cbind(1:maxit, rep(99, maxit) )
    error = FALSE

    
## * Find Tmat  ----
## Rotate F to PCA orientation P1 
  # we assume that P1 was orthogonally rotated to F
  findTmat <- function(F){
     W <- t(F) %*% F
     QDQ <- eigen(W, symmetric = TRUE)
     
     ## Q_FtoP1 is the matrix that takes F to P1 (unrotated components)
     Q <- QDQ$vectors
     # orient eigenvectors in positive direction
     Qsign <- sign(apply(Q, 2, sum))
     if(ncol(W) > 1) Qsign = diag(Qsign)
     Q <- Q %*% Qsign
     
     Q_FtoP1 <- Q
   
     
     ## Compute Tmat 
     # Tmat is the rotation matrix that takes P1 to F
     # Tmat <- as.matrix( solve(Q_FtoP1) ) since Q orthonormal
     Tmat <- as.matrix( t(Q_FtoP1) )
     
     list("Tmat" = Tmat,
          "Q_FtoP1" = Q_FtoP1)
  }#END findTmat   
  
## * Check if F legitimate PCA matrix ----
    checkF <- function(F){ 
      error = FALSE
      
      L1 <- eigen(F %*% t(F), symmetric = TRUE)$values
      sumL1 <- sum(L1)
      L1 <- L1[1:k]

      
      # eigenvalues must sum to p for legitimate F
      if (  (p-sumL1)/(p-k) > min(L1) ){
        error = TRUE
        warning("\n\nFATAL ERROR: F is not a legitimate PCA loadings matrix\n")
      }#END if (  (p-sumL1
      return(error)
    }#END checkF 
    
    
    
## * notConverged ----
##  Return values for non converged solution
    notConverged <- function(error, converged){
      return(
        list("R" = NULL,
             "Tmat" = NULL,
             "P1" = NULL,
             "Fhat" = NULL,
             "error" = error,
             "Lambda" = NULL,
             "iterHx" = iterHx,
             "converged" = converged,
             "Seed" = Seed)
      )#END return
    }#END notConverged   
    
  
#  * Barebones PCA function ----
  PCA <- function(R, k = NULL){
      if(is.null(k)) k <- ncol(R)
      VLV <- eigen(R, symmetric = TRUE)
      V <- VLV$vectors
      L <- VLV$values
  
      # Smooth if necessary
      L[L <.0] <- .01
  
     if( k > 1){
        P <-  V[, 1:k] %*% diag(L[1:k]^.5)
     }
     else{
        P <- as.matrix(V[, 1], drop=FALSE) * L[1]^.5
     }
      Psign <- sign(apply(P, 2, sum))
      if(k > 1) Psign = diag(Psign)
      P <- P %*%  Psign
      P
   }#END PCA   
  
# * InitP2Meth1----   
    # initialize P2 with a random bi-orthogonal matrix
    InitP2Meth1 <- function(P1, Seed){
      p <- nrow(P1)
      k <- ncol(P1)
      
      set.seed(Seed)
      P2 <- matrix( rnorm(p * (p - k)),
                    nrow = p,
                    ncol = (p-k))
    
      # orthogonalize P2
      P2 <- svd(P2)$u
 
      P2
    }#END InitP2Meth1
    
    
# * InitP2Meth2---- 
# initialize P2 with bi-orthogonal matrix 
# in null space of P1
    InitP2Meth2 <- function(P1, Seed){
      p <- nrow(P1)
      k <- ncol(P1)
      
      Ip <- diag(p)
      
      set.seed(Seed)
      P2 <- matrix( rnorm(p * (p - k)),
                    nrow = p,
                    ncol=(p-k))
      
      # Find P2 in Null space of P1
      P2 <- ( Ip - P1 %*% solve(t(P1) %*% P1) %*% t(P1)) %*% P2
      # orthogonalize P2
      P2 <- svd(P2)$u
      
      P2
    }#END PInitP2Meth2
    
    

# * PrjU ----
# Project on set of symmetric matrices with 
# unit diagonal  
  PrjU <- function(R){
    diag(R) <- 1
    R
  }  

  
# * PrjP ----
# Project on p x p matrices with P1 in the 
# leadings columns  
  PrjP <- function(R, P1, k, EigsSqrt){
    
    # If eigenvalues not supplied
    if(is.null(EigsSqrt)){
        P <- PCA(R)
        # replace leading PCA loadings with P1
        P[, 1:k] <- P1
    }
    
    # If eigenvalues are supplied
    if(!is.null(EigsSqrt)){
      VLV <- eigen(R)
      V2 <- VLV$vectors[,(k+1):p]
      # use fixed eigenvalues when creating P2
      P2 <- V2 %*% diag( EigsSqrt )
      P <- cbind( P1, P2)
    }
    P
}#END PrjP

# * testFit ----
  testFit <- function(R, Eigs){
    #if eigs not supplied then test fit by checking that
    # diag values of symmetric matrix  = 1
    if(is.null(Eigs)){
         #mean( abs(diag(R) - Ones_p) )
        fit <-  max( (diag(R) - Ones_p)^2 )
    }
    
    # if eigs are supplied then test eigs of symmetric matrix
    # and check that diag values = 1
    if(!is.null(Eigs)){
      L <- eigen(R, 
                 symmetric = T,
                 only.values = T)$values
      fit <- max( ( L - Eigs)^2 ) + max( (diag(R) - Ones_p)^2 )
      #fit <- sqrt(mean( ( L - Eigs)^2 ) )
    }
    
    fit
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # MAIN ----

  p <- nrow(F)
  k = ncol(F)
  Ones_p <- rep(1, p)
  
# * Step 1. Initialize P1----
  ## * Unrotate F to P1 
  fout <- findTmat(F)
  Tmat <- fout$Tmat
  Q_FtoP1 <- fout$Q_FtoP1
  P1 <- as.matrix( F %*% Q_FtoP1)
  
  # * Check if P1 is a legitimate PC loadings matrix ----
  Bad_P1 <- checkF(P1)

  # * If P1 ok then estimate R ----
  if(Bad_P1 == FALSE){
  
      set.seed(Seed)
      converged = TRUE
 
# * Step 1.5 Initialize P2----
      
      if(InitP2 == 1){
        # P2  = random semi orthogonal matrix
         P2 <- InitP2Meth1(P1, Seed)
      }
      if(InitP2 == 2){
        # P2 in left null space of P1
         P2 <- InitP2Meth2(P1, Seed)
      }   
      
  
# * Step 2. Initialize P----  
      P <- cbind(P1,P2)

# * Step 3. Execute MAP algorithm ----
      eps <- 999
      iter = 1
  
# R is a SSCP matrix
      R <- P %*% t(P)

      while(eps > epsMax){
    
      # Project on set of SSCP with unit diag
          R <- PrjU(R)
      
      # Project P on [P1 P2*]
          P <- PrjP(R, P1, k, EigsSqrt)
      
      #Compute implied SSCP and test fit
          R <- P %*% t(P)
      
      # ** Test fit ----
          eps <- testFit(R, Eigs)
      
      # Record iteration history
          iterHx[iter, 1] <- iter
          iterHx[iter, 2] <- eps
      
          if(PrintLevel == 2){
             cat("\n iter = ", iter, " eps = ",  eps )
          }  
      
          iter = iter + 1
          if(iter > maxit){
             converged = FALSE
             break
          } 
      
      }#END while   

      diag(R) <- 1
  
      P <- PCA(R)
      # if F was in rotated position then rotate P1 to F
      Fhat  <- P[, 1:k] %*% as.matrix(Tmat)
  
 # Only save iterations to convergence
      iterHx = iterHx[1:(iter-1), ]
 
      if(iter < maxit & PrintLevel >= 1){
         cat("\n\nProgram converged in ", iter, "iterations\n\n")
      } 
      if(iter > maxit & PrintLevel >= 1){
         cat("\n\nProgram failed to converge in ", iter, "iterations\n\n")
      } 
 
      Lambda <- eigen(R, symmetric = TRUE)$values
  } #End estimate R       
 

  # if P1 is not a legetimate loadings matrix return NULL
  if(Bad_P1 == TRUE){
    return(
    list("R" = NULL,
         "Tmat" = NULL,
         "P1" = NULL,
         "Fhat" = NULL,
         "error" = Bad_P1,
         "Lambda" = NULL,
         "iterHx" = NULL,
         "converged" = NULL,
         "Seed" = Seed)
    )
  }#END  error == TRUE 
    
 # Else return estimated results  
 ## Return ----
 if(Bad_P1 == FALSE){
    return(
    list("R" = R,
         "Tmat" = Tmat,
         "P1" = P[, 1:k],
         "Fhat" = Fhat,
         "error" = Bad_P1,
         "Lambda" = Lambda,
         "iterHx" = iterHx,
         "converged" = converged,
         "Seed" = Seed)
    )
  }#END  error == FALSE
  
}#END rPCA
  
  
  



