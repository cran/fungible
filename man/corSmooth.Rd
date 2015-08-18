\name{corSmooth}
\alias{corSmooth}
\title{Smooth a Non PD Correlation Matrix}
\description{
 A function for smoothing a non-positive definite correlation matrix by the method of Knol and Ten Berge (1989).
}
\usage{
corSmooth(R, eps = 1E8 * .Machine$double.eps)
}
\arguments{
 \item{R}{A non-positive definite correlation matrix.}
 \item{eps}{Small positive number to control the size of the smallest eigenvalue of the smoothed R matrix. Default = 1E8 * .Machine$double.eps}
}
\value{
\item{Rsmoothed}{A Smoothed (positive definite) correlation matrix.} 
}
\references{
Knol, D. and Ten Berge, J. (1989). Least-squares approximation of an improper correlation matrix by a proper one. \emph{Psychometrika, 54}, 53-61.
}
\author{Niels Waller}
\examples{
## choose eigenvalues such that R is NPD
l <- c(3.0749126,  0.9328397,  0.5523868,  0.4408609, -0.0010000)

## Generate NPD R
R <- genCorr(eigenval = l, seed = 123)
print(eigen(R)$values)

#> [1]  3.0749126  0.9328397  0.5523868  0.4408609 -0.0010000

## Smooth R
Rsm<-corSmooth(R, eps = 1E8 * .Machine$double.eps)
print(eigen(Rsm)$values)

#> [1] 3.074184e+00 9.326669e-01 5.523345e-01 4.408146e-01 2.219607e-08
}
\keyword{Statistics}% at least one, from doc/KEYWORDS
