\name{rMAP}
\alias{rMAP}
\title{Generate Correlation Matrices with Specified Eigenvalues}
\description{
  rMAP uses the method of alternating projections (MAP) to generate correlation matrices with specified eigenvalues.
}
\usage{
rMAP(eigenval, eps, maxits, Seed = NULL)
}
\arguments{
  \item{eigenval}{ A vector of eigenvalues that must sum to the order of the desired correlation matrix. A fatal error will occur if sum(eigenval) != length(eigenval).}
  \item{eps}{Convergence criterion. Default = 1e-12.}
  \item{maxits}{Maximm number of iterations of MAP.}
  \item{Seed}{Either a user supplied seed for the random number generator or  `NULL' for a function generated seed. Default Seed = `NULL'.}
}
\value{
\item{R}{A correlation matrix with the desired spectrum.}
\item{evals}{Eigenvalues of the returned matrix, R.}
\item{convergence}{(Logical) TRUE if MAP converged to a feasible solution, otherwise FALSE.}
}
\references{ 
Waller, N. G. (2016). Generating correlation matrices with specified eigenvalues using the method of alternating projections.
 }
\author{ Niels Waller}
\examples{

## Example
## Generate a correlation matrix with user-specified eigenvalues

R <- rMAP(c(2.5, 1, 1, .3, .2), Seed = 123)$R
print(R, 2)

#       [,1]    [,2]   [,3]    [,4]   [,5]
#[1,]  1.000  0.5355 -0.746 -0.0688 -0.545
#[2,]  0.535  1.0000 -0.671 -0.0016 -0.056
#[3,] -0.746 -0.6711  1.000  0.0608  0.298
#[4,] -0.069 -0.0016  0.061  1.0000  0.002
#[5,] -0.545 -0.0564  0.298  0.0020  1.000


eigen(R)$values
#[1] 2.5 1.0 1.0 0.3 0.2
}
\keyword{datagen}

