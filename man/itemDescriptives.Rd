\name{itemDescriptives}
\alias{itemDescriptives}
\title{Compute basic descriptives for binary-item analysis}
\description{
Compute basic descriptives for binary item analysis
}
\usage{
itemDescriptives(X, digits=3)
}
\arguments{
 \item{X}{a matrix of binary (0/1) item responses.}
 \item{digits}{number of digits to print.}
}
\value{ 
\item{alpha}{Coefficient alpha for the total scale.}
\item{means}{item means.}
\item{standard deviations}{item standard deviations.}
\item{pt. biserial correlations}{corrected item-total point biserial correlations.}
\item{biserial correlations}{corrected item-total point biserial correlations.}
\item{corrected.alpha}{corrected (leave item out) alpha coefficients.}
}
\author{Niels Waller}
\examples{
	## Example 1: generating binary data to match
	## an existing binary data matrix
	##
	## Generate correlated scores using factor 
	## analysis model
	## X <- Z *L' + U*D 
	## Z is a vector of factor scores
	## L is a factor loading matrix
	## U is a matrix of unique factor scores
	## D is a scaling matrix for U

	Nsubj <- 2000
	L <- matrix( rep(.707,5), nrow = 5, ncol = 1)
	Z <-as.matrix(rnorm(Nsubj))
	U <-matrix(rnorm(Nsubj * 5),nrow = Nsubj, ncol = 5)
	tmp <-  sqrt(1 - L^2) 
	D<-matrix(0, 5, 5)
	diag(D) <- tmp
	X <- Z \%*\% t(L) + U\%*\%D

	cat("\nCorrelation of continuous scores\n")
	print(round(cor(X),3))

	thresholds <- c(.2,.3,.4,.5,.6)

	Binary<-matrix(0,Nsubj,5)
	for(i in 1:5){
	  Binary[X[,i]<=thresholds[i],i]<-1
	}   

	cat("\nCorrelation of Binary scores\n")
	print(round(cor(Binary),3))

	## Now use 'bigen' to generate binary data matrix with 
	## same correlations as in Binary

	z <- bigen(data = Binary, n = 5000)

	cat("\n\nnames in returned object\n")
	print(names(z))

	cat("\nCorrelation of Simulated binary scores\n")
	print(round( cor(z$data), 3))


	cat("Observed thresholds of simulated data:\n")
	cat( apply(z$data, 2, mean) )
	
	itemDescriptives(z$data)
}

\keyword{statistics}% at least one, from doc/KEYWORDS
