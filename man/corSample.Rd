\name{corSample}
\alias{corSample}
\title{Sample Correlation Matrices from a Population Correlation Matrix}
\description{
Sample correlation (covariance) matrices from a population correlation matrix (see Browne, 1968; Kshirsagar, 1959)
}
\usage{
corSample(R, n)
}
\arguments{
 \item{R}{A population correlation matrix.}
 \item{n}{Sample correlation (covariance) matrices will be generated assuming a sample size of n.}
}
\value{
\item{cor.sample}{Sample correlation matrix.} 
\item{cov.sample}{Sample covariance matrix.}
}
\references{
Browne, M. (1968). A comparison of factor analytic techniques. \emph{Psychometrika, 33(3)}, 267-334.
 
Kshirsagar, A. (1959). Bartlett decomposition and Wishart distribution. \emph{The Annals of Mathematical Statistics, 30(1)}, 239-241.
}
\author{Niels Waller}
\examples{
R <- matrix(c(1, .5, .5, 1), 2, 2)
# generate a sample correlation from pop R with n = 25
out <- corSample(R, n = 25)
out$cor.sample
out$cov.sample
}
\keyword{datagen}% at least one, from doc/KEYWORDS
