\name{kurt}
\alias{kurt}
\title{Calculate Univariate Kurtosis for a Vector or Matrix}
\description{
Calculate univariate kurtosis for a vector or matrix (algorithm G2 in Joanes & Gill, 1998).
}
\usage{
kurt(x)
}
\arguments{
  \item{x}{Either a vector or matrix of numeric values.}
}
\value{
  \item{Kurtosis for each column in x.}{}
}
\references{
Joanes, D. N. & Gill, C. A. (1998). Comparing measures of sample skewness and kurtosis. \emph{The Statistician, 47}, 183-189.
}
\author{Niels Waller}
\seealso{ 
  \code{\link{skew}}
}
\examples{
x <- matrix(rnorm(1000), 100, 10)
print(kurt(x))
}
\keyword{Statistics}% at least one, from doc/KEYWORDS


