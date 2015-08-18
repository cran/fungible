\name{monte1}
\alias{monte1}
\title{Simulate Multivariate Non-normal Data by Vale & Maurelli (1983) Method}
\description{
Function for simulating multivariate nonnormal data by the methods described by Fleishman (1978) and Vale & Maurelli (1983).
}
\usage{
monte1(seed, nvar, nsub, cormat, skewvec, kurtvec)
}
\arguments{
\item{seed}{An integer to be used as the random number seed.}
\item{nvar}{Number of variables to simulate.}
\item{nsub}{Number of simulated subjects (response vectors).}
\item{cormat}{The desired correlation matrix.} 
\item{skewvec}{A vector of indicator skewness values.}
\item{kurtvec}{A vector of indicator kurtosis values.}
}
\value{
\item{data}{The simulated data.}
\item{call}{The call.}
\item{nsub}{Number of subjects.}
\item{nvar}{Number of variables.}
\item{cormat}{The desired correlation matrix.}
\item{skewvec}{The desired indicator skewness values.}
\item{kurtvec}{The desired indicator kurtosis values.}
\item{seed}{The random number seed.}
}
\references{
Fleishman, A. I (1978). A method for simulating non-normal distributions. \emph{Psychometrika, 43}, 521-532.

Vale, D. C., & Maurelli, V. A. (1983). Simulating multivariate nonnormal distributions. \emph{Psychometrika, 48}, 465-471.

} 
\author{Niels Waller}
\seealso{\code{\link{monte}}, \code{\link{summary.monte}}, \code{\link{summary.monte1}}
}
\examples{

## Generate dimensional data for 4 variables. 
## All correlations = .60; all variable
## skewness = 1.75; 
## all variable kurtosis = 3.75
 
cormat <- matrix(.60,4,4)
diag(cormat) <- 1

nontaxon.dat <- monte1(seed = 123, nsub = 100000, nvar = 4, skewvec = rep(1.75, 4),
               kurtvec = rep(3.75, 4), cormat = cormat)
 
print(cor(nontaxon.dat$data), digits = 3)
print(apply(nontaxon.dat$data, 2, skew), digits = 3)
print(apply(nontaxon.dat$data, 2, kurt), digits = 3)               
}
\keyword{datagen}% at least one, from doc/KEYWORDS
