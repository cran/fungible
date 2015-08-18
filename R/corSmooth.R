
##--------------------------------------#
##  October 31, 2012
##  updated August 12, 2015
##  Niels Waller
##
##  Smooth a non-positive definate R matrix by the method of Knol and Ten Berge
##--------------------------------------#

corSmooth <- function(R, eps = 1E8 * .Machine$double.eps){
  
     ULU <- eigen(R)
     L <- ULU$values
     L[L<=0]<- eps
     
     U <- ULU$vectors
     L<-ncol(R)*L/sum(L)
     R<-U %*% diag(L) %*% t(U)
     D<-diag(1/sqrt(diag(R)))
     D%*%R%*%D
}

