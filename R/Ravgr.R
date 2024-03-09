# March 09, 2024
# 
#' Generate a random R matrix with an average rij
#'
#' Ravgr(Rseed, NVar = NULL, u = NULL, rdist = "U", alpha = 4, beta = 2, SEED = NULL)
#' 
#' @param Rseed  (matrix or scalar) This argument can take one of two alternative inputs. 
#' The first input is an \eqn{n \times n} \strong{R} matrix with a known, 
#' average rij.  The second type of input is a scalar \eqn{\bar{r}_{ij}}. 
#' @param NVar (integer)  If \code{Rseed} is a scalar then the user must specify
#' \code{NVar}, the number of variables in the desired \strong{R} matrix. Default(\code{NVar = NULL}).
#' @param u (scalar). A scalar \eqn{\in [0,1]}. Higher values of u will produce
#' \strong{R} matrices with more variable off-diagonal elements. 
#' @param rdist (character). A character that controls the variance of the off 
#' diagonal elements of the generated \strong{R}. If \code{u = NULL} and \code{rdist="U"} then the \strong{R} 
#' matrices are uniformly  sampled from the space of all \eqn{n\times n} \strong{R} matrices with a 
#' fixed average rij.  If \code{u = NULL} and \code{rdist = "B"} then the \strong{R} matrices are 
#' selected as a function of the \code{alpha} and \code{beta} arguments of a 
#' Beta distribution. Default \code{rdist= "U"}. See Waller (2024) for details.
#' @param alpha (numeric) The shape1 parameter of a beta distribution.
#' @param beta (numeric)  The shape2 parameter of a beta distribution.
#' @param SEED (numeric)  The initial seed for the random number generator. 
#' If SEED is not supplied then the program will generate (and return) a randomly
#' generated seed.
#'
#' @return 
#' \itemize{
#'  \item \strong{R} A random \strong{R} matrix with a known, average off-diagonal 
#'  element rij.
#'  \item \strong{Rseed} The input \strong{R} matrix or scalar with the
#'  desired average rij.
#'  \item \strong{u} A user-specified or random number \eqn{\in [0,1]}.
#'  \item \strong{s} Scaling factor for hollow matrix \code{H}.
#'  \item \strong{H}  A hollow matrix used to create a fungible \strong{R} matrix.
#'  \item \strong{alpha}  First argument of the beta distribution. If 
#'  \code{rdist= "U"} then \code{alpha = NULL}.
#'  \item \strong{beta}  Second argument of the beta distribution. 
#'  If \code{rdist= "U"} then \code{beta = NULL}.
#'  \item  \strong{SEED} The initial value for the random number generator.
#'}  
#' 
#' 
#' @references Waller, N. G. (2024). Generating correlation matrices with a 
#' user-defined average correlation. Manuscript under review. 
#'
#' @author Niels G. Waller
#'
#' @examples
#'  # Example 1
#'   R <- matrix(.35, 6, 6)
#'   diag(R) <- 1
#'   
#'   Rout <- Ravgr(Rseed = R, 
#'                rdist = "U", SEED = 123)$R
#'                
#'   Rout |> round(3)            
#'   mean( Rout[upper.tri(Rout, diag = FALSE)] )
#'   
#'   # Example 2 
#'   Rout <- Ravgr(Rseed = .35, NVar = 6, 
#'                rdist = "U", SEED = 123)$R
#'                
#'   Rout |> round(3)            
#'   mean( Rout[upper.tri(Rout, diag = FALSE)] )   
#'   
#'   # Example 3
#'   # Generate an R matrix with a larger var(rij)
#'   Rout <- Ravgr(Rseed = .35,
#'                NVar = 6, 
#'                rdist = "B",
#'                alpha = 7,
#'                beta = 2)$R
#'                
#'   Rout |> round(3)            
#'   mean( Rout[upper.tri(Rout, diag = FALSE)] )
#'   
#'   # Example 4: Demonstrate the function of u
#'   sdR <- function(R){
#'     sd(R[lower.tri(R, diag = FALSE)])
#'   }
#'   
#'   Rout <- Ravgr(Rseed = .35,
#'                NVar = 6, 
#'                u = 0,
#'                SEED = 123)
#'   sdR(Rout$R)  
#'   
#'   Rout <- Ravgr(Rseed = .35,
#'                NVar = 6, 
#'                u = .5,
#'                SEED = 123)
#'   sdR(Rout$R)   
#'   
#'   Rout <- Ravgr(Rseed = .35,
#'                NVar = 6, 
#'                u = 1,
#'                SEED = 123)
#'   sdR(Rout$R)          
#'    
#' @export
#' 
#' 
Ravgr <- function(Rseed, NVar = NULL, u = NULL, rdist = "U", 
                  alpha = 4, beta = 2, SEED = NULL){
  # if a = 1 and b = 1 then u ~ U(0,1)
  
  R = Rseed
  
  # Check input arguments 
     if(!is.matrix(Rseed) & !is.null(NVar)) {
        R <- matrix(Rseed, NVar, NVar)
        diag(R) <- 1
        Rseed <- R
     }
  
    if(!is.matrix(Rseed) & is.null(NVar)){
      stop("\n\nFATAL ERROR: Invalid entry for NVar")
    }
  
   if(is.matrix(Rseed)){
     NVar <- ncol(Rseed)
   } #END of input checks
  
  if(!is.null(u)){
    if(u < 0 | u > 1) stop("\n\nFATAL ERROR: Invalid entry for u")
  }
  
  ## generate random seed if not supplied
  if(is.null(SEED)) SEED <- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) 
  set.seed(SEED)
  
  # Define vecp.inv
  vecp.inv <- function(x, NVar){
     M <- matrix(0, NVar, NVar)
     M[upper.tri(M, diag = FALSE)] <- x
     M = M + t(M)
     M
  }# END vecp.inv
  

  # generate vector of standardized deviation scores 
  Delta <- scale(rnorm( .5* NVar * (NVar - 1), 0, 1) )
  
  # H is a hollow matrix
  H <- vecp.inv(Delta, NVar)

  # Find s s.t. Rk is on the surface of the elliptope
  f <- function(s){
    (eigen(R + s * H, symmetric = TRUE, only.values = TRUE)$val[NVar])^2
  }
    
   s <- optimize(f, c(0,100),
             maximum = FALSE,
             tol = 1E-8)$minimum
  
  
   if(is.null(u)){
    # generate u: a random number in [0,1] 
     if(rdist =="U") u <- runif(1, 0, 1)
     if(rdist =="B") u <- rbeta(1, alpha, beta)
   
     if(rdist == "U"){
      alpha = NULL
      beta  = NULL
     }
   }##END if(is.null(u)) 
  
  # Rk is a fungible R matrix with avg rij 
  Rk <- R + u*s*H
  
  
  # Return Values
  invisible(
    list(
       R = Rk,  # a fungible R matrix with avg rij
       Rseed = Rseed,
       u = u,
       s = s,
       H = H,
       alpha = alpha,
       beta = beta,
       SEED = SEED
       ))
}

