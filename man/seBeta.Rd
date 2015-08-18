\name{seBeta}
\alias{seBeta}
\title{Standard Errors and CIs for Standardized Regression Coefficients}
\description{
Computes Normal Theory and ADF Standard Errors and CIs for Standardized Regression Coefficients
}
\usage{
seBeta(X, y, cov.x = NULL, cov.xy = NULL, var.y = NULL, Nobs = NULL, 
alpha = 0.05, ADF = TRUE, digits = 3)		
 }
\arguments{
 \item{X}{Matrix of predictor scores.}  
 \item{y}{Vector of criterion scores.}
 \item{cov.x}{Covariance or correlation matrix of predictors.}
 \item{cov.xy}{Vector of covariances or correlations between predictors and criterion.}
 \item{var.y}{Criterion variance.}
 \item{Nobs}{Number of observations.}
 \item{alpha}{Desired Type I error rate; default = .05.}
 \item{ADF}{Logical (TRUE/FALSE) to select ADF confidence intervals - requires raw X and raw y; default = TRUE. If ADF = FALSE then the function will compute Normal Theory confidence intervals.}
 \item{digits}{Number of significant digits to print; default = 3.}
}
\value{
\item{cov.Beta}{Normal theory or ADF covariance matrix of standardized regression coefficients.} 
\item{se.Beta}{standard errors for standardized regression coefficients.} 
\item{alpha}{desired Type-I error rate.}  
\item{CI.Beta}{Normal theory or ADF (1-alpha)\% confidence intervals for standardized regression coefficients.}
\item{estimator}{estimator = "ADF" or "Normal".} 
}

\references{
Jones, J. A, and Waller, N. G. (2015). The Normal-Theory and Asymptotic 
Distribution-Free (ADF) covariance matrix of standardized regression coefficients: Theoretical extensions and finite sample behavior. \emph{Psychometrika, 80}, 365-378.
}
\author{Jeff Jones and Niels Waller}
\examples{
library(MASS)
R <- matrix(.5, 3, 3)
diag(R) <- 1
X <- mvrnorm(n = 200, mu  =rep(0, 3), Sigma = R, empirical = TRUE)
Beta <- c(.2, .3, .4)
y <- X \%*\% Beta + .64 * scale(rnorm(200))
seBeta(X, y, Nobs = 200, alpha = .05, ADF = TRUE)
}

\keyword{Statistics}% at least one, from doc/KEYWORDS
