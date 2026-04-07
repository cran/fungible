# April 6, 2026
#
#' Generate a random R matrix with given MAD from a user-defined, population R
#'
#' Rmad(Rpop, MAD = .10, tries = 1E5, SEED = NULL)
#' 
#' @param Rpop  (matrix) A model-implied, population correlation matrix. 
#' @param MAD (scalar) The desired mean absolute deviation value (i.e., 
#' \eqn{mad(Rpop - R) = MAD} where \eqn{mad} computes the mean absolute deviation 
#' for the nonredundant, non-diagonal values of its argument).
#' @param tries (scalar) The maximum number of tries for the optimization algorithm. 
#' Default(\eqn{tries = 1E5}).
#' @param SEED (numeric)  The starting seed for the random number generator. 
#' If SEED is not supplied then the program will generate (and return) a randomly
#' generated seed.
#'
#' @return 
#' \itemize{
#'  \item \strong{R} A random \strong{R} matrix with a known MAD from \strong{Rpop}.
#'  \item \strong{Rpop} The model-implied, population correlation matrix.
#'  \item \strong{Delta} The vector of model-errors.
#'  \item \strong{H}  A hollow matrix used to create a fungible \strong{R} matrix.
#'  \item \strong{iter} Iteration number at convergence (or max iteration number).
#'  \item \strong{MAD} The target mean absolute deviation.
#'  \item \strong{converged}. A Boolean that describes the convergence status of 
#'  the optimization algorithm.
#'  \item  \strong{SEED} The initial value for the random number generator.
#'}  
#' 
#'
#' @author Niels Waller
#'
#' @examples
#'  # Example 1
#'   Rpop <- matrix(.35, 6, 6)
#'   diag(Rpop) <- 1
#'   
#'   out <- Rmad(Rpop, MAD = .038, SEED = 123)
#'                
#'   out$R |> round(3)            
#'  
#' # E = mattix of model error values
#'   E = Rpop - out$R 
#  # Compute MAD
#'   mean(abs(E[lower.tri(E, diag = FALSE)]))
#'    
#' @export
#' 
#' 
Rmad <- function(Rpop, MAD = .03, tries= 1E5, SEED = NULL){
  
  Rpop  #user-supplied input matrix
  NVar = ncol(Rpop)
  
  ## Generate random seed if not supplied
  if(is.null(SEED)) SEED <- as.integer((as.double(Sys.time())*10000+Sys.getpid()) %% 2^31) 
  set.seed(SEED)
  
  # Define vecp.inv
  #  See Waller, N. G. (2024). A simple and fast algorithm for 
  #  generating correlation matrices with a known average correlation coefficient. 
  #  The American Statistician. 
  vecp.inv <- function(x, NVar){
     M <- matrix(0, NVar, NVar)
     M[upper.tri(M, diag = FALSE)] <- x
     M = M + t(M)
     M
  }# END vecp.inv
  
 # getMAD:  
 # Find a solution that yields the desired MAD
 # For justification of this method, see Chapter 5 in:
 # Devroye, L. (1986). Non-Uniform Random Variate Generation. 
 # New York, Berlin, Heidelberg Tokyo: Springer-Verlag. 
  
 # We want uniformly sampled points that lie on the intersection
 # of an elliptope (the space of R) and the surface of a cross polytope 
 # (i.e., a high-dimensional diamond like shape)
  getMAD <- function(Rpop){
    
    # Select points uniformly on the surface of the cross polytope
    m <- NVar * (NVar - 1)/ 2
    U <- runif(m-1, 0, 1)
    
    V <- sort(U, decreasing = FALSE)
    V <- c(0, V, 1)
    
    #Scale the vector to a fixed L1 norm
    X <- m * MAD * diff(V)
    
    signs <- sample(c(1, -1), size = m, replace = TRUE)
    
    # Delta is sampled from the desired density
    Delta <- X * signs
  
    # Insert Delta into H, a hollow matrix
    H <- vecp.inv(Delta, NVar)
     
   # Rk (when PSD) is a random R matrix with the desired MAD
    Rk <- Rpop + H
    
    list( Rk = Rk,
          Delta = Delta,
          H = H)
  }#END getMAD
  
  converged = FALSE
  iter = 0
  #cat("\nThinking very hard . . . ")

  for(i in 1:tries){
    iter = iter + 1
    SEED <- SEED + 1
    set.seed(SEED)
    # search for a feasible solution
    out = getMAD(Rpop)
    Rk = out$Rk
    # Check that Rk is PSD
    if( eigen(Rk, only.values = TRUE, symmetric = TRUE)$values[NVar] >= 0 ){
      #cat("\nSolution found in", iter, "iterations")
      converged = TRUE
      break
    }#END if  
  }#END for (i in 1:tries)
  
  if(iter >= tries){
    cat("\n*** ERROR: No feasible solution found ***")
  }
  
  # Return Values
  invisible(
    list(
       R = Rk,  # a fungible R matrix with the desired MAD
       Rpop = Rpop,
       Delta = out$Delta,
       H = out$H,
       iter = iter,
       MAD = MAD,
       converged = converged,
       SEED = SEED
       ))
}



