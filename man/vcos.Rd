\name{vcos}
\alias{vcos}
\title{Compute the Cosine Between Two Vectors}
\description{
Compute the cosine between two vectors. 
}
\usage{
vcos(x,y)
}
\arguments{
 \item{x}{A p x 1 vector.}
 \item{y}{A p x 1 vector.}
}
\value{
 \item{Cosine between x and y}{} 
}
\examples{
x <- rnorm(5)
y <- rnorm(5)
vcos(x, y)
}
\keyword{Statistics}% at least one, from doc/KEYWORDS
