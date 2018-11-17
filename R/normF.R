normF <- function(X){
## return the Frobenius norm of the input matrix 
  ## trace fnc 
  ## TR <- function(X) sum(diag(X))
  ## sqrt(TR(t(X) %*% X))
  ## Frobenius norm
  sqrt(sum(X^2))
}  