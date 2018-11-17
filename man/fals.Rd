\name{fals}
\alias{fals}
\title{Unweighted least squares factor analysis}
\description{
Unweighted least squares factor analysis
}
\usage{
fals(R,  
     nfactors, 
     TreatHeywood = TRUE) 
}
\arguments{
 \item{R}{Input correlation matrix.}
 \item{nfactors}{Number of factors to extract.} 
 \item{TreatHeywood}{If TreatHeywood = TRUE then a penalized least squares function is used to bound the commonality estimates below 1.0.  Default(TreatHeywood = TRUE).}
}
\value{
 \item{loadings}{Unrotated factor loadings. If a Heywood case is present in the initial solution then the model is re-estimated via non-iterated principal axes with max(rij^2) as fixed communaility (h2) estimates.}
 \item{h2}{Vector of final commonality estimates.}
 \item{uniqueness}{Vector of factor uniquenesses, i.e. (1 - h2).}
 \item{Heywood}{(logical) TRUE if a Heywood case was produced in the LS solution.}
 \item{TreatHeywood}{(logical) Value of the TreatHeywood argument.}
 \item{converged}{(logical) TRUE if all values of the gradient are sufficiently close to  zero.}
 \item{MaxAbsGrad}{The maximum absolute value of the gradient at the solution.
  }
}
\author{Niels Waller}
\examples{

Rbig <- fungible::rcor(120)
out1 <- fals(R = Rbig, 
             nfactors = 2,
             TreatHeywood = TRUE)
}
\keyword{Statistics}% at least one, from doc/KEYWORDS
