
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
    QLQ <- eigen(Su)
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
  
  evals <- eigen(R)$values
  list(R=R, evals=evals, convergence=convergence)
} ## END rMAP





###################################################################
##              
###################################################################

# ev<-c(2,.75,.25)
# rMAP(ev)
