\name{rmsd}
\alias{rmsd}
\title{
Root Mean Squared Deviation of (A - B) 
}
\description{
Calculates the root mean squared deviation of
matrices A and B.  If these matrices are symmetric (Symmetric = TRUE)
then the calculation is based on the upper triangles of each 
matrix. When the matrices are symmetric, the diagonal of each 
matrix can be included or excluded from the calculation 
(IncludeDiag = FALSE)
}
\usage{
rmsd(A, B, Symmetric = TRUE, IncludeDiag = FALSE)
}

\arguments{
  \item{A}{A possibly non square matrix.}
  \item{B}{A matrix of the same dimensions as  matrix A.}
  \item{Symmetric}{Logical indicating whether A and B are symmetric matrices. (Default: Symmetric = TRUE)}
  \item{IncludeDiag}{Logical indicating whether to include the diagonals
  in the calculation. (Default: IncludeDiag = FALSE).}
}

\value{
Returns the root mean squared deviation of (A - B). 
}
\author{
Niels Waller
}

\examples{
A <- matrix(rnorm(9), nrow = 3)
B <- matrix(rnorm(9), nrow = 3)

( rmsd(A, B, Symmetric = FALSE, IncludeDiag = TRUE) )
}
\keyword{stats}

