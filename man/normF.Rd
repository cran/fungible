\name{normF}
\alias{normF}
\title{Compute the Frobenius norm of a matrix}
\description{A function to compute the Frobenius norm of a matrix}
\usage{ 
normF(X)
}
\arguments{
\item{X}{A matrix.}
}
\value{
\item{}{The Frobenius norm of X.}
}
\author{Niels Waller}
\examples{
data(BadRLG)
out <- smoothLG(R = BadRLG, Penalty = 50000)
cat("\nGradient at solution:", out$gr,"\n")
cat("\nNearest Correlation Matrix\n")
print( round(out$RLG,8) )
cat("\nFrobenius norm of (NPD - PSD) matrix\n")
print(normF(BadRLG - out$RLG ))
}
\keyword{statistics}% at least one, from doc/KEYWORDS

