BiFAD <- function(R, B = NULL, nGroup = NULL, 
                  factorMethod = "minres", 
                  rotation="oblimin",  salient = .25, 
                  maxitFA = 5000, 
                  maxitRotate = 5000,
                  gamma = 0){
  
  ###############################################################
  ## AUTHOR: Niels Waller
  ## August 18, 2017
  ##  requires  psych:   package for initial factor extraction  
  ##            GPArotation:  for rotation options
  ##  
  ## Arguments:
  ##  R:  Input correlation matrix
  ##  
  ##  B:  bifactor target matrix. If B=NULL the program will 
  ##        create an empirically defined target matrix. 
  ##
  ##  nGroup: Number of group factors in bifactor solution.
  ##
  ##  factorMethod: factor extraction method. Options include:  
  ##        minres (minimum residual), ml (maximum likelihood), 
  ##        pa (principal axis), gls (generalized least squares). 
  ##  
  ##  rotation:  factor rotation method. Current options include: 
  ##        oblimin, geominQ, quartimin, promax.  
  ##  
  ##  salient:  Threshold value for creating an empirical target 
  ##        matrix.
  ##  
  ##  maxitFA:  Maximum iterations for the factor extraction 
  ##        method.
  ##  
  ##  maxitRotate: Maximum iterations for the gradient pursuit 
  ##        rotation algorithm.
  ##  
  ##  gamma: Optional tuning parameter for oblimin rotation.
  ##
  ##  Value:
  ##  
  ##    B: User defined or empirically generated target matrix.
  ##  
  ##    BstarSL:  Direct S-L solution.
  ##
  ##    BstarFR:  Direct full rank bifactor solution.
  ##  
  ##    rmsrSL:   Root mean squared residual of (B - BstarSL). 
  ##
  ##    rmsrFR:   Root mean squared residual of (B - BstarFR).  
  ################################################################  
  
  ## Compute orthogonal Procrustes rotation
  Procrustes <-function(M1, M2){
    tM1M2 <- t(M1) %*% M2
    svdtM1M2 <- svd(tM1M2)
    P <- svdtM1M2$u
    Q <- svdtM1M2$v
    T <- Q %*% t(P)
    ## Orthogonally rotate M2 to M1
    M2 %*% T
  }   
  
  ## FATAL ERROR
  if(is.null(B) & is.null(nGroup)){
    stop("\n\n FATAL ERROR: Either B or nGroup must be specified") 
  }  
  
  ## Compute unrotated factor solution
  if(is.null(nGroup))  nGroup <- ncol(B) - 1 
  
  F <- psych::fa(r = R, 
                 nfactors = nGroup, 
                 fm = factorMethod, 
                 rotate = "none")$loadings[]
  
  ## Append column of zeros to create rank deficient matrix
  L0 <- cbind(F,0)
  
  
  ## Create 0/1 Target matrix if B not supplied  
  Bflag <- 1 # set to 1 if Target matrix (B) is known
  if(is.null(B))  Bflag <- 0
  
  if(!Bflag & (!is.null(nGroup))){ 
    
    # Start rotation from random position
    TmatRandom <- qr.Q(qr(matrix(rnorm(nGroup * nGroup ), 
                                 nGroup )))
    switch(rotation, 
           "oblimin" = {
             gpa.out <- GPArotation::oblimin(F, Tmat=TmatRandom, 
                                             gam = gamma, 
                                             maxit=maxitRotate)
           }, "geominQ" = {
             gpa.out <- GPArotation::geominQ(F, 
                                             Tmat = TmatRandom,
                                             maxit=maxitRotate)
           }, "quartimin" = {
             gpa.out <- GPArotation::quartimin(F, 
                                               Tmat = TmatRandom, 
                                               maxit=maxitRotate)
           }, "promax" = {
             gpa.out <- promax(F)
           })
    
    B <- gpa.out$loadings[]
    
    ## Record signs of loadings
    signB <- sign(B)
    
    ## Convert Target matrix into signed 0/1 matrix
    B <- signB * matrix(as.numeric(abs(B) >= salient), 
                       nrow(R), ncol = nGroup)
    
    ## Append ones vector for general factor
    B <- cbind(1,B)
    
  } ## END Create 0/1 Target matrix
  
  
  ## Rotate rank deficient loading matrix to Target
  BstarSL <- Procrustes(B, L0)
  
  ## Compute non-hierarchical bifactor solution
  F2    <- psych::fa(r = R, 
                     nfactors = nGroup+1, 
                     fm=factorMethod, 
                     rotate="none")$loadings[]
  
  ## Rotate full rank loading matrix to Target
  BstarFR <- Procrustes(B, F2)
  
  ## Compute root mean squared residual of Target matrix (B)
  ## and best fitting SL and FR solutions (Bstar)
  rmsrSL <- rmsrFR <- NA
  if(Bflag == 1){
    rmsrSL = sqrt(mean((B - BstarSL)^2))
    rmsrFR = sqrt(mean((B - BstarFR)^2))
  }   
  
  list(B = B,           
       BstarSL = BstarSL,
       BstarFR = BstarFR,
       rmsrSL = rmsrSL,
       rmsrFR = rmsrFR
  )
} ## END BiFAD

