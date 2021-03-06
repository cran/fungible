% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genPhi.R
\name{genPhi}
\alias{genPhi}
\title{Create a random Phi matrix with maximum factor correlation}
\usage{
genPhi(NFac, EigenValPower = 6, MaxAbsPhi = 0.5)
}
\arguments{
\item{NFac}{Number of factors.}

\item{EigenValPower}{(Scalar > 1) A scalar than controls the positive
skewness of the distribution of eigenvalues of Phi.}

\item{MaxAbsPhi}{(Scaler in [0,1]) The maximum off diagonal of Phi (the
factor correlation matrix).}
}
\value{
A factor correlation matrix. Note that the returned matrix is not guaranteed 
to be positive definite. However, a PD check is performed in simFA so that simFA always 
produces a PD Phi matrix.
}
\description{
Create a random Phi matrix with maximum factor correlation.
}
\examples{

NFac <- 5
par(mfrow=c(2,2))
  for(i in 1:4){
     R <- genPhi(NFac, 
               EigenValPower = 6, 
               MaxAbsPhi = 0.5)
               
    L <- eigen(R)$values
    plot(1:NFac, L, 
        type="b",
        ylab = "Eigenvalues of Phi",
        xlab = "Dimensions",
        ylim=c(0,L[1]+.5))
  }

}
\author{
Niels Waller
}
\keyword{stats}
