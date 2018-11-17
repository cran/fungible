
fals   <- function(R,           ## Correlation matrix
                  nfactors,
                  TreatHeywood = TRUE){    ## Number of factors to extract
  
 

  # vers: November 16, 2018
  # Author: Niels Waller
  #
  # Purpose:       Compute an unweighted least squares factor analysis
  #
  # Args:
  #    R:          Correlation matrix
  #    nfactors:   Number of factors to extract
  #   TreatHeywood: If True a penality is applied
  #          to the function to avoid heywood cases
  #
  # Output:
  #    loadings:   Matrix of unrotated factor loadings
  #       if a Heywood case was present in the initial solution then the 
  #       model is re-estimated via non-iterated principal axses with max(rij^2) 
  #       as fixed communaility (h2) estimates.
  #    h2:         Vector of communality estimates 
  #    uniqueness: Vector of factor uniquenesses (1 - h2)
  #    Heywood: (logical)  TRUE if a Heywood case was produced
  #                in the LS solution.
  #   converged: (logical) TRUE is gradient is zero
  #   MaxAbsGrad:  The maximum absolute value of the gradient at the solution
  
  NVar <- nrow(R)
  Nparam <- NVar * nfactors
  
 
  ## ---- Calculate OLS discrepancy function ----
  LS.f <- function(Fmat){
    Fmat <- matrix(Fmat, 
                   nrow  = NVar, 
                   ncol  = nfactors, 
                   byrow = FALSE)
    
    
    Sigma <- tcrossprod(x = Fmat, 
                        y = Fmat)
    
    diag(Sigma) <- 1
  
    
    #LS fit function
    sum(  (Sigma - R)^2 ) 
    
  } #END Calculate OLS discrepancy function 
  
  
  ## ---- Calculate Penalized OLS discrepancy function ----
  PenalizedLS.f <- function(Fmat){
    Fmat <- matrix(Fmat, 
                   nrow  = NVar, 
                   ncol  = nfactors, 
                   byrow = FALSE)
    
    
    Sigma <- tcrossprod(x = Fmat, 
                        y = Fmat)
    
    diag(Sigma) <- 1
    
    penalty = 1
    if( max(apply(Fmat^2,1,sum)) >= .99) penalty = 5E4
    #penalized LS fit function
    penalty * sum(  (Sigma - R)^2 ) 
    
  } #END Calculate OLS discrepancy function 
  
  
  # ---- Calculate gradient of OLS discrepancy function---- 
  grdFALS <-function(Fmat){
    Fmat <- matrix(Fmat, 
                   nrow  = NVar, 
                   ncol  = nfactors, 
                   byrow = FALSE)
    
    Sigma <- tcrossprod(x = Fmat, 
                        y = Fmat)
    
    Usq <- diag(1 - diag(Sigma))
    
    ## November 13, 2018
    ## variances must be non negative
    Usq[Usq < 0 ] <- 0
    
    # See Mulaik 2nd edition EQ 8.43
    4 * tcrossprod(x = Fmat, y = Fmat) %*% 
      Fmat - 4 * R %*% 
      Fmat  + 4 * 
      crossprod(x = Usq, y = Fmat)
  } #END Calculate gradient
  
  
  
  # ---- Generate start values via non-iterated PA solution ---- 
    Rstart <- R
    diag(Rstart) <- 0
  
    # Use max rij^2 for each col
    maxrij <- apply(abs(Rstart), 2, max)
    diag(Rstart) <- maxrij^2
  
  
    if(NVar < 40){
      ULU <- eigen(Rstart)
      
      if(nfactors == 1){
        Fstart <- as.matrix(ULU$vectors[, 1] * sqrt(ULU$values[1]))
      }
      
      if(nfactors > 1) {
        ULU <- eigen(Rstart)   
        #abs(eigval) added for possible npd
        Fstart <- ULU$vectors[, 1:nfactors] %*% 
        diag(sqrt(abs(ULU$values[1:nfactors])))
      }
      
      
      # Treat Heywood cases in Start Values
      initialH2 <- apply(Fstart^2,1,sum)
      if(max(initialH2) > 1){
        s <- sqrt(rep(.99, NVar)/h2)
        Fstart <-diag(s) %*% Fstart 
      }
    
      
    }## End if (NVar < 40)
  
    if(NVar >= 40) {
      # calculate the k largest eigenvalues and associated eigenvectors 
      ULU <- RSpectra::eigs(Rstart, 
                  k     = nfactors, 
                  which = "LM", 
                  lower = TRUE)
    
      if(nfactors == 1){
        Fstart <- as.matrix(ULU$vectors[, 1] * sqrt(ULU$values[1]))
      }
      if(nfactors > 1) {
        Fstart <- ULU$vectors[,1:nfactors] %*% 
          diag(sqrt(abs(ULU$values[1:nfactors])))
      }
      
      # Treat Heywood cases in Start Values
      initialH2 <- apply(Fstart^2,1,sum)
      if(max(initialH2) > 1){
        s <- sqrt(rep(.80, NVar)/h2)
        Fstart <-diag(s) %*% Fstart 
      }
      
     
    }## end if
  
    # ---- minimize Non Penalized OLS function----  
    if(TreatHeywood == FALSE){  
      nParam <- NVar * nfactors  
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      # ---- Calculate least squares solution via "BFGS" ----
      output <- optim(par     = Fstart, 
                      fn      = LS.f, 
                      gr      = grdFALS, 
                      method  =  "L-BFGS-B",
                      lower =  rep(-1,  nParam),
                      upper =  rep( 1,  nParam),
                      control = list(maxit=1E5))
      
      # fl = unrotated factor loadings
      fl <- matrix(output$par, 
                   nrow  = NVar,
                   ncol  = nfactors,
                   byrow = FALSE) 
      
    } #END if(TreatHeywood == FALSE)
    
    
    
  # ---- minimize Penalized OLS function----  
  if(TreatHeywood == TRUE){  
     nParam <- NVar * nfactors  
     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
     # Use BFGS in initial estimation  
     # ---- Calculate least squares solution via "BFGS" ----
     output <- optim(par     = Fstart, 
                    fn      = PenalizedLS.f, 
                    gr      = grdFALS, 
                    method  =  "L-BFGS-B",
                    lower =  rep(-1,  nParam),
                    upper =  rep( 1,  nParam),
                    control = list(maxit=1E5))
                                   
     # fl = unrotated factor loadings
     fl <- matrix(output$par, 
               nrow  = NVar,
               ncol  = nfactors,
               byrow = FALSE) 
  
  } #END if(TreatHeywood = TRUE)
 

    #####################################
    # check new solution for Heywood case
    h2 <- apply(fl^2, 1, sum)
    maxh2 <- max(h2) 
    Heywood <- FALSE
    if(maxh2 > 1) Heywood <-TRUE
    
  
    uniqueness <- 1 - h2
  
    # ---- Compute max absolute gradient element at solution ----
    grdFALS <- grdFALS(fl)
    MaxAbsGrad <- max( abs(grdFALS ) )
  
  list(loadings   = fl, 
       h2         = h2, 
       uniqueness = uniqueness,
       Heywood    = Heywood,
       converged = as.logical(1 - output$convergence),
       MaxAbsGrad = MaxAbsGrad,
       grdFALS = grdFALS)

} #end function fals

#########################################################


