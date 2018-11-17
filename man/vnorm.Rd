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
 the scaled (i.e., unit length) input vector
}
\author{Niels Waller}
\examples{
 x <- rnorm(5)
 v <- vnorm(x)
 print(v)
}
\keyword{statistics } 
