#' Compute the Frobenius norm of a matrix
#' 
#' A function to compute the Frobenius norm of a matrix
#' 
#' 
#' @param X A matrix.
#' @return \item{}{The Frobenius norm of X.}
#' @author Niels Waller
#' @keywords statistics
#' @export
#' @examples
#' 
#' data(BadRLG)
#' out <- smoothLG(R = BadRLG, Penalty = 50000)
#' cat("\nGradient at solution:", out$gr,"\n")
#' cat("\nNearest Correlation Matrix\n")
#' print( round(out$RLG,8) )
#' cat("\nFrobenius norm of (NPD - PSD) matrix\n")
#' print(normF(BadRLG - out$RLG ))
#' 
normF <- function(X){
## return the Frobenius norm of the input matrix 
  ## trace fnc 
  ## TR <- function(X) sum(diag(X))
  ## sqrt(TR(t(X) %*% X))
  ## Frobenius norm
  sqrt(sum(X^2))
}  
