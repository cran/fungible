\name{skew}
\alias{skew}
\title{Calculate Univariate Skewness for a Vector or Matrix}
\description{
Calculate univariate skewness for vector or matrix (algorithm G1 in Joanes & Gill, 1998).
}
\usage{
skew(x)
}
\arguments{
  \item{x}{Either a vector or matrix of numeric values.}
}
\value{
 \item{Skewness for each column in x.}{}
}
\references{
Joanes, D. N. & Gill, C. A. (1998). Comparing measures of sample skewness and kurtosis. \emph{The Statistician, 47}, 183-189.
}
\author{Niels Waller}
\seealso{\code{\link{kurt}}
}
\examples{
x <- matrix(rnorm(1000), 100, 10)
skew(x)
}
\keyword{Statistics}% at least one, from doc/KEYWORDS

