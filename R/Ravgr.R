# April 28, 2023
# 
#' Generate a random R matrix with an average rij
#'
#' Ravgr(R, rdist = "U", alpha = 4, beta = 2)
#' 
#' @param R  (matrix) An n x n correlation matrix with a known, average rij.
#' @param rdist (character). A character that controls the variance of the off 
#' diagonal elements of the generated R. If \code{rdist="U"} then the R 
#' matrices are uniformly  sampled from the space of all n x n R matrices with a 
#' fixed average rij.  If \code{rdist = "B"} then the R matrices are 
#' selected as a function of the \code{alpha} and \code{beta} arguments of a 
#' Beta distribution. Default \code{rdist= "U"}. See Waller (2023) for details.
#' @param alpha (numeric) The shape1 parameter of a beta distribution.
#' @param beta (numeric)  The shape2 parameter of a beta distribution.
#' 
#'
#' @return A random R matrix with a known, average off-diagonal element.
#' 
#' @references Waller, N. G. (2023). Generating correlation matrices with a 
#' user-defined average correlation. Manuscript under review. 
#'
#' @author Niels G. Waller
#'
#' @examples
#'
#'   R <- matrix(.35, 6, 6)
#'   diag(R) <- 1
#'   
#'   Rout <- Ravgr(R, 
#'                rdist = "U")
#'                
#'   Rout |> round(3)            
#'   mean( Rout[upper.tri(Rout, diag = FALSE)] )
#'   
#'   
#'   # Generate an R matrix with a larger var(rij)
#'   Rout <- Ravgr(R, 
#'                rdist = "B",
#'                alpha = 7,
#'                beta = 2)
#'                
#'   Rout |> round(3)            
#'   mean( Rout[upper.tri(Rout, diag = FALSE)] )
#'    
#' @export
#' 
#' 
Ravgr <- function(R, rdist = "U", alpha = 4, beta = 2){
  # if a = 1 and b = 1 then u ~ U(0,1)
  
  NVar <- ncol(R)

  vecp.inv <- function(x, NVar){
     M <- matrix(0, NVar, NVar)
     M[upper.tri(M, diag = FALSE)] <- x
     M = M + t(M)
     M
  }# END vecp.inv
  

  # generate vector of deviation scores 
  Delta.star <- rnorm( .5* NVar * (NVar - 1), 0, 1) 
  # center scores
  Delta <- Delta.star - mean(Delta.star)
  
  # H is a hollow matrix
  H <- vecp.inv(Delta, NVar)

  # Find s s.t. Rk is on the surface of the elliptope
  f <- function(s){
    (eigen(R + s * H, symmetric = TRUE, only.values = TRUE)$val[NVar])^2
  }
    
   s <- optimize(f, c(0,100),
             maximum = FALSE,
             tol = 1E-8)$minimum
  
  # generate a random number in [0,1] 
   if(rdist =="U") u <- runif(1, 0, 1)
   if(rdist =="B") u <- rbeta(1, alpha, beta)
  
  Rk <- R + u*s*H
  Rk
}

