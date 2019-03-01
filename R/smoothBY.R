#' Smooth an NPD R matrix to PD using the Bentler Yuan 2011 method
#' 
#' Smooth a NPD correlation matrix to PD using the Bentler and Yuan method.
#' 
#' 
#' @param R Indefinite Matrix.
#' @param const const is a user-defined parameter that is defined as k in
#' Bentler and Yuan (2011). If 0 < const < 1, then const is treated as a fixed
#' value. If const = 1 then the program will attempt to find the highest value
#' of const such that R is positive (semi) definite.
#' @param eps If const = 1 then the program will iteratively reduce const by
#' eps until either (a) the program converges or (b) const < = 0.
#' @return \item{RBY}{smoothed correlation matrix.} \item{constant}{The final
#' value of const.} \item{convergence}{(Logical) a value of TRUE indicates that
#' the function converged.} \item{outStatus}{Convergence state for Rcsdp::csdp
#' function. \cr \cr 0: \cr \cr Success. Problem solved to full accuracy \cr
#' \cr 1: \cr \cr Success. Problem is primal infeasible \cr \cr 2: \cr \cr
#' Success. Problem is dual infeasible \cr \cr 3: \cr \cr Partial Success.
#' Solution found but full accuracy was not achieved \cr \cr 4: \cr \cr
#' Failure. Maximum number of iterations reached \cr \cr 5: \cr \cr Failure.
#' Stuck at edge of primal feasibility \cr \cr 6: \cr \cr Failure. Stuch at
#' edge of dual infeasibility \cr \cr 7: \cr \cr Failure. Lack of progress \cr
#' \cr 8:\cr \cr Failure. X or Z (or Newton system O) is singular \cr \cr 9:
#' \cr \cr Failure. Detected NaN or Inf values} \cr \item{glb}{Greatest lower
#' bound reliability estimates.} \item{eps}{Default value (eps = 1E-03) or
#' user-supplied value of eps.}
#' @author Code modified from that reported in Debelak, R. & Tran, U. S.
#' (2011).
#' @references Bentler, P. M. & Yuan, K. H.  (2011).  Positive definiteness via
#' off-diagonal scaling of a symmetric indefinite matrix.  \emph{Psychometrika,
#' 76}(1), 119--123.
#' 
#' Debelak, R. & Tran, U. S. (2013). Principal component analysis of smoothed
#' tetrachoric correlation matrices as a measure of dimensionality.
#' \emph{Educational and Psychological Measurement, 73}(1), 63--77.
#' @keywords statistics
#' @import Rcsdp
#' @export
#' @examples
#' 
#' data(BadRBY)
#' 
#' out<-smoothBY(R = BadRBY, const = .98)
#' cat("\nSmoothed Correlation Matrix\n")
#' print( round(out$RBY,8) )
#' cat("\nEigenvalues of smoothed matrix\n")
#' print( eigen(out$RBY)$val  )
#' 
smoothBY <- function(R, const = .98, eps = .001){
  ## smoothBY is modified from code originally published in
  ## Debelak, R. & Tran, U. S. (2013). Principal component analysis of 
  ## smoothed tetrachoric correlation matrices as a measure of dimensionality. 
  ## Educational and Psychological Measurement, 73(1), 63--77. 
  
  ## A function for smoothing an improper correlation matrix to a 
  ## psd matrix based on theory of Bentler and Yuan 2011.
  
  ## Arguments
  ## R        indefinite "pseudo" correlation matrix
  ## const    const is defined as k in Bentler and Yuan. If const is a numeric 
  ##          such that 0 < constant < 1, then const is treated as a fixed value.
  ##          If const = 1 then the program will attempt to find the highest 
  ##          value of const such that R is positive (semi) definite.
  ##
  ## Value
  ## RBY     smoothed correlation matrix
  ## constant The final value of const, 
  ## convergence (Logical) TRUE if smoothBY converged to feasible solution, 
  ##             False otherwise. 
  ##
  ## Update: March 11, 2018 (added const argument)
  ## Update: March 29, 2018 fixed code to replicate Bentler & Yuan 
  ## Update: April 01, 2108 argument eps added
  
   if (!requireNamespace("Rcsdp")) {
        stop("Rcsdp must be installed to run smoothBY")
    }

   if(const < 0 || const > 1)  stop("const must be between 0 and 1")
  
   
   Nvar <- ncol(R)
   
   # initialize glb
   glb <- rep(999, Nvar)
   
   # Upbounds =  av in Debelack and Tran 2013
   UpBounds <- as.vector(rep(2.5,Nvar))
   LoBounds <- rep(0, Nvar)
  
  

    opt = rep(1, Nvar)
    
   ##  C:	 A list defining the block diagonal cost matrix C.
    C <- list(diag(Nvar) - R, -UpBounds, LoBounds)

   ##  A: 	A list of length m containing block diagonal constraint 
   ##      matrices A_i. Each constraint matrix A_i is specified by a 
   ##      list of blocks as explained in the Details section.    
    A <- vector("list", Nvar)

        for (i in 1:Nvar) {
   ## b:	A numeric vector of length m containing the right hand side of 
   ##     the constraints.      
        b <- rep(0, Nvar)
        b[i] <- 1
        A[[i]] <- list(diag(b), -b, b)
    }
    
   ## K: Describes the domain of each block of 
   ##    the sdp problem. It is a list with the following elements:
   ## type:
   ##   A character vector with entries "s" or "l" indicating the type 
   ##   of each block. If the jth entry is "s", then the jth block is a 
   ##   positive semidefinite matrix. otherwise, it is a vector with 
   ##   non-negative entries.
   ## size:
   ##   A vector of integers indicating the dimension of each block.  
    K <- list(type = c("s", "l", "l"), size = rep(Nvar, 3))
    
    out <- try(Rcsdp::csdp(C, A, opt, K, 
                           control = Rcsdp::csdp.control(printlevel = 0, maxiter=200)),
                           silent = FALSE)
    
    # cat("\n\n")
    # print(out$status)
    
      if (class(out) == "try-error" ) {
        warning("Rcsdp (try-error) convergence failure: ", out$status)
        return(list(RBY = NULL, 
                    constant = const, 
                    convergence = FALSE,
                    outStatus = out$status,
                    glb = glb,
                    eps = eps))
      }
      
      if(class(out) !="try-error"){
        if( out$status!= 0){
          warning("Rcsdp Convergence failure: ", out$status)
          return(list(RBY = NULL,
                      constant = const,
                      convergence = FALSE,
                      outStatus = out$status,
                      glb = glb,
                      eps = eps))
        }
        RcsdpConvergence <- out$status
      }
  
  
  
    ## RcsdpConvergence == 0 indicates problem solved to full accuracy
  if(RcsdpConvergence == 0){
    ## y is the optimal dual solution. A vector of the same length as b
  estimatedCommunalities <- glb <- out$y

   # As reported by the original authors:  
   # This code calculates the communalities of each 
   # variable after applying a minimum trace factor analysis as 
   # values of the vector ew. The used functions require the 
   # definition of upper limits for these communalities. In 
   # Debelack and Tran, these limits were defined as 2.5. The number of 
   # items is denoted by Nvar. In a second step, the correlations 
   # are rescaled using a scaling constant. 
  
  # Find optimal value of const
  minEig <- -99
  convergence = TRUE
  
  ## const supplied by user
  if( is.numeric(const) & const >= 0  & const < 1) {
      RBY <- R
      for (i in 1:Nvar) {
        if (estimatedCommunalities[i] >= 1.0) { 
          for (j in 1:Nvar) {
            if(i!=j) {
              RBY [i,j] <- RBY [j,i] <- RBY [i,j] * sqrt(const/estimatedCommunalities[i]) 
            } # end if (i!=i)
          } # end for b in 1:Nvar
        } # end if estimated
      } # end for a in 1:Nvar
      
      
      # Convergence check
      # is the smallest eigenvalue less than 0
      minEig <- eigen(RBY)$val[Nvar]
      if( (minEig < 0) || (max(abs(RBY)) > 1)){
        warning("Convergence failure: Lmin < 0  or |rij| > 1")
        return(list(RBY = NULL, 
                    constant = const, 
                    convergence = FALSE,
                    outStatus = out$status,
                    glb = glb,
                    eps = eps))
      }
      else
      return(list(RBY = RBY, 
                  constant = const, 
                  convergence = convergence, 
                  outStatus = out$status,
                  glb = glb,
                  eps = eps))
  } ## End if is numeric(const)
  
  if(const == 1) {
  while(minEig < 0){
    RBY <- R
    for (i in 1:Nvar) {
       if (estimatedCommunalities[i] >= 1.0) { 
        for (j in 1:Nvar) {
          if(i!=j) {
          RBY [i,j] <- RBY [j,i] <- RBY [i,j] * sqrt(const/estimatedCommunalities[i]) 
        } # end if (i!=i)
       } # end for b in 1:Nvar
     } # end if estimated
    } # end for a in 1:Nvar
    
    # is the smallest eigenvalue less than 0
    minEig <- eigen(RBY)$val[Nvar]
    if( minEig < 0 ){
       const <- const - eps
    }   
    
    # Check if: 0 < const < 1 
    if(const < 0 ){
      warning("Convergence failure: const < 0")
      return(list(RBY = NULL, constant = NULL, convergence = FALSE, eps = eps))
    }
    
  } ## END of: while minEig < 0
    
    # Check if all |RBY_{ij}| < = 1
    if(max(abs(RBY)) > 1){
      warning("Convergence failure: |rij| > 1")
      list(RBY = NULL, constant = NULL, convergence = FALSE, eps = eps)
    } 
    
    return(list(RBY = RBY, 
                constant = const, 
                convergence = convergence,
                outStatus = out$status,
                glb = glb,
                eps = eps))  
  
   } ## END if const == 1
  
 
  } ## End of: if RcsdpConvergence == 0 (converged)
    
} # END smoothBY 
###################################################################

