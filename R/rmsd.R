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