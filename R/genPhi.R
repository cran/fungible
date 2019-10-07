# Create random Phi with max phi_ij



#' Create a random Phi matrix with maximum factor correlation
#' 
#' Create a random Phi matrix with maximum factor correlation.
#' 
#' 
#' @param NFac Number of factors.
#' @param EigenValPower (Scalar > 1) A scalar than controls the positive
#' skewness of the distribution of eigenvalues of Phi.
#' @param MaxAbsPhi (Scaler in [0,1]) The maximum off diagonal of Phi (the
#' factor correlation matrix).
#' @return A factor correlation matrix. Note that the returned matrix is not guaranteed 
#' to be positive definite. However, a PD check is performed in simFA so that simFA always 
#' produces a PD Phi matrix.
#' @author Niels Waller
#' @keywords stats
#' @examples
#' 
#' NFac <- 5
#' par(mfrow=c(2,2))
#'   for(i in 1:4){
#'      R <- genPhi(NFac, 
#'                EigenValPower = 6, 
#'                MaxAbsPhi = 0.5)
#'                
#'     L <- eigen(R)$values
#'     plot(1:NFac, L, 
#'         type="b",
#'         ylab = "Eigenvalues of Phi",
#'         xlab = "Dimensions",
#'         ylim=c(0,L[1]+.5))
#'   }
#' 
#' @export genPhi
genPhi <- function(NFac, EigenValPower = 6, MaxAbsPhi =.5){
  # generate random positive numbers that sum to NFac
  eigs <- sort(abs(rnorm(NFac))^EigenValPower, decreasing=TRUE)
  eigs <- eigs * NFac/sum(eigs)
  L <- diag(eigs)
  
  # generate random orthogonal matrices using QR decomposition
  M <- matrix(rnorm(NFac*NFac), nrow=NFac, ncol=NFac)
  Q <- qr.Q(qr(M))
  Cov <- Q %*% L %*% t(Q)
  # scale to R matrix
  R <- cov2cor(Cov)
  # make hollow matrix
  R <- R - diag(NFac)
  maxR <- max(abs(R))
  S <- diag(NFac) * sqrt(MaxAbsPhi)/sqrt(maxR)
  R <- S %*% R %*% S
  diag(R) <- 1
  Phi <- round(R,2)
  Phi
}






