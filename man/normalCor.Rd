\name{normalCor}
\alias{normalCor}
\title{Compute Normal-Theory Covariances for Correlations}
\description{
 Compute normal-theory covariances for correlations
}
\usage{
normalCor(R, Nobs)
}
\arguments{
 \item{R}{a p x p matrix of correlations.}
 \item{Nobs}{Number of observations.}
}
\value{
A normal-theory covariance matrix of correlations. 
}
\references{
Nel, D.G. (1985). A matrix derivation of the asymptotic covariance matrix of sample correlation coefficients. \emph{Linear algebra and its applications, 67}, 137--145.
}
\author{Jeff Jones and Niels Waller}
\seealso{\code{\link{adfCor}}}
\examples{
	data(Harman23.cor)
	normalCor(Harman23.cor$cov, Nobs = 305)
}

\keyword{Statistics}% at least one, from doc/KEYWORDS
