\name{rcor}
\alias{rcor}
\title{Generate Random PSD Correlation Matrices}
\description{
Generate random PSD correlation matrices.
}
\usage{
rcor(Nvar)
}
\arguments{
  \item{Nvar}{An integer that determines the order of the random correlation matrix.}
}
\value{
  \item{A random correlation matrix.}{}
}
\details{
 rcor generates random PSD correlation matrices by (1) generating Nvar squared random normal deviates, (2) scaling the deviates to sum to Nvar, and then (3) placing the scaled values into a diagonal matrix L. Next, (4) an Nvar x Nvar orthogonal matrix, Q, is created by performing a QR decomposition of a matrix, M, that contains random normal deviates.  (5) A PSD covariance matrix, C, is created from Q L Q^T and then (6) scaled to a correlation metric.
 }
\author{Niels Waller}
\seealso{ 
  \code{\link{genCorr}}
}
\examples{
R <- rcor(4)
print( R )
}
\keyword{Statistics}% at least one, from doc/KEYWORDS


