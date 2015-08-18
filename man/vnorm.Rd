\name{vnorm}
\alias{vnorm}
\title{Norm a Vector to Unit Length}
\description{Norm a vector to unit length.}
\usage{
vnorm(x)
}
\arguments{
 \item{x}{An n by 1 vector.}
}
\value{
 \item{x scaled to unit length.}{}
}
\author{Niels Waller}
\examples{
x <- rnorm(5)
v <- vnorm(x)
}
\keyword{Statistics}

