###################################################################
# rcor
# function to generate random PSD R matrices
# August 10, 2016
# Author: Niels Waller
#
# Arguments:
#     Nvar   An integer that determines the order of the correlation matrix.
# Output
#            An Nvar x Nvar PSD correlation matrix.

rcor <- function(Nvar){
  
  # generate random positive numbers that sum to Nvar
  eigs <- rnorm(Nvar)^2
  eigs <- eigs * Nvar/sum(eigs)
  L <- diag(rev(sort(eigs)))
  
  # generate random orthogonal matrices using QR decomposition
  M <- matrix(rnorm(Nvar*Nvar), nrow=Nvar, ncol=Nvar)
  Q <- qr.Q(qr(M))
  Cov <- Q %*% L %*% t(Q)
  
  # scale to R matrix
  cov2cor(Cov)
}