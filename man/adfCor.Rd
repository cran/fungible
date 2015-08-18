\name{adfCor}
\alias{adfCor}
\title{Asymptotic Distribution-Free Covariance Matrix of Correlations}
\description{
  Function for computing an asymptotic distribution-free covariance matrix of 
  correlations.
}
\usage{
adfCor(X, y = NULL)
}
\arguments{
 \item{X}{Data matrix.}
 \item{y}{Optional vector of criterion scores.}
}
\value{
\item{adfCorMat}{Asymptotic distribution-free estimate of the covariance matrix of correlations.}
}
\references{
Browne, M. W. (1984). Asymptotically distribution-free methods for the 
	analysis of covariance structures. \emph{British Journal of Mathematical and Statistical 
	Psychology, 37}, 62--83.
	
Steiger, J. H. and Hakstian, A. R. (1982). The asymptotic distribution of elements
	of a correlation matrix: Theory and application. \emph{British Journal of Mathematical 
	and Statistical Psychology, 35}, 208--215.
}
\author{Jeff Jones and Niels Waller}
\examples{
## Generate non-normal data using monte1
set.seed(123)
## we will simulate data for 1000 subjects
N <- 1000

## R = the desired population correlation matrix among predictors
R <- matrix(c(1, .5, .5, 1), 2, 2)

## Consider a regression model with coefficient of determination (Rsq): 
Rsq <- .50

## and vector of standardized regression coefficients
Beta <- sqrt(Rsq/t(sqrt(c(.5, .5))) \%*\% R \%*\% sqrt(c(.5, .5))) * sqrt(c(.5, .5))

## generate non-normal data for the predictors (X)
## x1 has expected skew = 1 and kurtosis = 3
## x2 has expected skew = 2 and kurtosis = 5
X <- monte1(seed = 123, nvar = 2, nsub = N, cormat = R, skewvec = c(1, 2), 
            kurtvec = c(3, 5))$data
            
## generate criterion scores            
y <- X \%*\% Beta + sqrt(1-Rsq)*rnorm(N)

## Create ADF Covariance Matrix of Correlations
adfCor(X, y)

#>             12           13           23
#> 12 0.0012078454 0.0005331086 0.0004821594
#> 13 0.0005331086 0.0004980130 0.0002712080
#> 23 0.0004821594 0.0002712080 0.0005415301

}

\keyword{Statistics}

