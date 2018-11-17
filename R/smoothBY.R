
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

