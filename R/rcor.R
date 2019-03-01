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



#' Generate Random PSD Correlation Matrices
#' 
#' Generate random PSD correlation matrices.
#' 
#' rcor generates random PSD correlation matrices by (1) generating Nvar
#' squared random normal deviates, (2) scaling the deviates to sum to Nvar, and
#' then (3) placing the scaled values into a diagonal matrix L. Next, (4) an
#' Nvar x Nvar orthogonal matrix, Q, is created by performing a QR
#' decomposition of a matrix, M, that contains random normal deviates.  (5) A
#' PSD covariance matrix, C, is created from Q L Q^T and then (6) scaled to a
#' correlation metric.
#' 
#' @param Nvar An integer that determines the order of the random correlation
#' matrix.
#' @return \item{A random correlation matrix.}{}
#' @author Niels Waller
#' @seealso \code{\link{genCorr}}
#' @keywords Statistics
#' @export
#' @examples
#' 
#' R <- rcor(4)
#' print( R )
#' 
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
