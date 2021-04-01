#' Generate Correlation Matrices with Specified Eigenvalues
#' 
#' rMAP uses the method of alternating projections (MAP) to generate
#' correlation matrices with specified eigenvalues.
#' 
#' 
#' @param eigenval A vector of eigenvalues that must sum to the order of the
#' desired correlation matrix. A fatal error will occur if sum(eigenval) !=
#' length(eigenval).
#' @param eps Convergence criterion. Default = 1e-12.
#' @param maxits Maximm number of iterations of MAP.
#' @param Seed Either a user supplied seed for the random number generator or
#' `NULL' for a function generated seed. Default Seed = `NULL'.
#' @return \item{R}{A correlation matrix with the desired spectrum.}
#' \item{evals}{Eigenvalues of the returned matrix, R.}
#' \item{convergence}{(Logical) TRUE if MAP converged to a feasible solution,
#' otherwise FALSE.}
#' @author Niels Waller
#' @references Waller, N. G. (2016). Generating correlation matrices with
#' specified eigenvalues using the method of alternating projections.
#' @keywords datagen
#' @export
#' @examples
#' 
#' 
#' ## Example
#' ## Generate a correlation matrix with user-specified eigenvalues
#' 
#' R <- rMAP(c(2.5, 1, 1, .3, .2), Seed = 123)$R
#' print(R, 2)
#' 
#' #       [,1]    [,2]   [,3]    [,4]   [,5]
#' #[1,]  1.000  0.5355 -0.746 -0.0688 -0.545
#' #[2,]  0.535  1.0000 -0.671 -0.0016 -0.056
#' #[3,] -0.746 -0.6711  1.000  0.0608  0.298
#' #[4,] -0.069 -0.0016  0.061  1.0000  0.002
#' #[5,] -0.545 -0.0564  0.298  0.0020  1.000
#' 
#' 
#' eigen(R)$values
#' #[1] 2.5 1.0 1.0 0.3 0.2
#' 
rMAP <- function (eigenval, eps = 1e-12, maxits = 5000, Seed = NULL) 
{
  ###################################################################
  ## rMAP: Generate a correlation matrix (R) with a fixed set of 
  ##       eigenvalues by the method of alternating projections 
  ##
  ## Niel Waller 
  ## September 24, 2016
  ## October 27, 2016
  ##
  ## Input:
  ## eigenval    Desired spectrum of R
  ## eps         convergence criterion. Default = 1e-12
  ## maxits       Maximum number of iterations of the MAP
  ## Seed        Optional seed for generating initial 
  ##             orthogonal matrix
  ## Output:
  ## R           A correlation matrix with the desired spectrum
  ## evals       Eigenvalues of R at solution
  ## convergence (Logical) TRUE if MAP converged to a feasible solution, 
  ##             otherwise FALSE
  ###################################################################
  
  ## February 11, 2021:  To speed up code, substitue the following 
  ## when executing the quadratic forms (LambdaStar should be 
  ## a vector not a diagonal matrix)
  ## S.LambdaStar <- crossprod(t(Q) * sqrt(LambdaStar))
  
  Nvar <- length(eigenval)
  if (!isTRUE(all.equal(sum(eigenval), Nvar))) 
    stop("Sum of eigenvalues not equal to Number of Variables\n")
  
  ## generate random seed if not supplied
  if(is.null(Seed)) Seed<- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) 
  set.seed(Seed)
  
  ## generate a random orthogonal matrix
  M <- matrix(rnorm(Nvar * Nvar), nrow = Nvar, ncol = Nvar)
  Q <- qr.Q(qr(M))
  
  ## create a PSD covariance matrix with desired spectrum
  
   LambdaStar <- diag(eigenval)
   S.LambdaStar <- Q %*% LambdaStar %*% t(Q)
   
  
  # enforce symmetry
  S.LambdaStar <- .5 * (S.LambdaStar + t(S.LambdaStar))
  
  delta <- 1
  iter <- 0
  while(delta >= eps){

    ## Project onto symmetric matrix with unit diagonals
    Su <- S.LambdaStar
    diag(Su) <- 1
    
    #update eigenvectors
    QLQ <- eigen(Su, symmetric = TRUE)
    Q <- QLQ$vectors
    Lambda <- QLQ$values
    
    ## test convergence
    delta <- sqrt(sum((eigenval-Lambda)^2))
    iter <- iter + 1
    if(iter > maxits) {
      warning("\nFailed to converge after maxits iterations")
      delta = 0
    }  
    ## Project onto symmetric matrix with given spectrum
    S.LambdaStar <- Q %*% LambdaStar %*% t(Q)
    # enforce symmetry
    S.LambdaStar <- .5 * (S.LambdaStar + t(S.LambdaStar))
  }
  
  # at convergence
  R <- S.LambdaStar
  diag(R) <- 1
  
  # check solution feasibility 
  convergence <- FALSE
  if(max(abs(R[upper.tri(R)])) <= 1 & iter < maxits) convergence<-TRUE
  
  evals <- eigen(R, symmetric = TRUE)$values
  list(R=R, evals=evals, convergence=convergence)
} ## END rMAP





###################################################################
##              
###################################################################

# ev<-c(2,.75,.25)
# rMAP(ev)
