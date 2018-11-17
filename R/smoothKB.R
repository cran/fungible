


smoothKB <- function(R, eps = 1E8 * .Machine$double.eps){
  ##--------------------------------------#
  ##  October 31, 2012
  ##  updated August 12, 2015
  ##  updated  Sept 25 2015
  ##  Niels Waller
  ##
  ##  Smooth a non-positive definite R matrix by the method of 
  ##  Knol and Berger 1991 Eq (27)
  ##--------------------------------------#
  
     KDK <- eigen(R)
     Dplus <- D <- KDK$values
     Dplus[D <= 0]<- eps
     Dplus <- diag(Dplus)
     K <- KDK$vectors
    
     Rtmp <- K %*% Dplus %*% t(K)
     invSqrtDiagRtmp <- diag(1/sqrt(diag(Rtmp)))
     RKB <- invSqrtDiagRtmp %*% Rtmp %*% invSqrtDiagRtmp
      
    colnames(RKB) <- rownames(RKB) <- colnames(R)

	 list(RKB = RKB, eps = eps)
}

