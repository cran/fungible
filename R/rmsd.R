##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## FNC: rmsd
## AU: Niels Waller
## DATE: January 16, 2018
##
## Description: Calculates the root mean squared deviation of
## two matrices.  If the matrices are symmetric (Symmetric = TRUE)
## then the calculation is based on the upper triangles of each 
## matrix. When the matrices are symmetric, the diagonal of each 
## matrix can be included or excluded from the calculation 
## (IncludeDiag = FALSE)
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



#' Root Mean Squared Deviation of (A - B)
#' 
#' Calculates the root mean squared deviation of matrices A and B.  If these
#' matrices are symmetric (Symmetric = TRUE) then the calculation is based on
#' the upper triangles of each matrix. When the matrices are symmetric, the
#' diagonal of each matrix can be included or excluded from the calculation
#' (IncludeDiag = FALSE)
#' 
#' 
#' @param A A possibly non square matrix.
#' @param B A matrix of the same dimensions as matrix A.
#' @param Symmetric Logical indicating whether A and B are symmetric matrices.
#' (Default: Symmetric = TRUE)
#' @param IncludeDiag Logical indicating whether to include the diagonals in
#' the calculation. (Default: IncludeDiag = FALSE).
#' @return Returns the root mean squared deviation of (A - B).
#' @author Niels Waller
#' @keywords stats
#' @export
#' @examples
#' 
#' A <- matrix(rnorm(9), nrow = 3)
#' B <- matrix(rnorm(9), nrow = 3)
#' 
#' ( rmsd(A, B, Symmetric = FALSE, IncludeDiag = TRUE) )
#' @export

rmsd <- function(A,
                 B,
                 Symmetric = TRUE,
                 IncludeDiag = FALSE) {
  if (Symmetric) {
    sqrt(mean((A[lower.tri(A, diag = IncludeDiag)] -
                 B[lower.tri(B, diag = IncludeDiag)])^2))
  }
  else{
    sqrt(mean((A - B)^2))
  }
}
