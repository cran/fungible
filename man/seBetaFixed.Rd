% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seBetaFixed.R
\name{seBetaFixed}
\alias{seBetaFixed}
\title{Covariance Matrix and Standard Errors for Standardized Regression
Coefficients for Fixed Predictors}
\usage{
seBetaFixed(
  X = NULL,
  y = NULL,
  cov.x = NULL,
  cov.xy = NULL,
  var.y = NULL,
  var.error = NULL,
  Nobs = NULL
)
}
\arguments{
\item{X}{Matrix of predictor scores.}

\item{y}{Vector of criterion scores.}

\item{cov.x}{Covariance or correlation matrix of predictors.}

\item{cov.xy}{Vector of covariances or correlations between predictors and
criterion.}

\item{var.y}{Criterion variance.}

\item{var.error}{Optional argument to supply the error variance: var(y -
yhat).}

\item{Nobs}{Number of observations.}
}
\value{
\item{cov.Beta}{Normal theory covariance matrix of standardized
regression coefficients for fixed predictors.} \item{se.Beta}{Standard
errors for standardized regression coefficients for fixed predictors.}
}
\description{
Computes Normal Theory Covariance Matrix and Standard Errors for
Standardized Regression Coefficients for Fixed Predictors
}
\examples{

## We will generate some data and pretend that the Predictors are being held fixed

library(MASS)
R <- matrix(.5, 3, 3); diag(R) <- 1
Beta <- c(.2, .3, .4)

rm(list = ".Random.seed", envir = globalenv()); set.seed(123)
X <- mvrnorm(n = 200, mu = rep(0, 3), Sigma = R, empirical = TRUE)
y <- X \%*\% Beta + .64*scale(rnorm(200))

seBetaFixed(X, y)

# $covBeta
#              b1           b2           b3
# b1  0.003275127 -0.001235665 -0.001274303
# b2 -0.001235665  0.003037100 -0.001491736
# b3 -0.001274303 -0.001491736  0.002830157
# 
# $seBeta
#         b1         b2         b3 
# 0.05722872 0.05510989 0.05319922

## you can also supply covariances instead of raw data

seBetaFixed(cov.x = cov(X), cov.xy = cov(X, y), var.y = var(y), Nobs = 200)

# $covBeta
#              b1           b2           b3
# b1  0.003275127 -0.001235665 -0.001274303
# b2 -0.001235665  0.003037100 -0.001491736
# b3 -0.001274303 -0.001491736  0.002830157
# 
# $seBeta
#         b1         b2         b3 
# 0.05722872 0.05510989 0.05319922


}
\references{
Yuan, K. & Chan, W. (2011). Biases and standard errors of
standardized regression coefficients. \emph{Psychometrika, 76(4)}, 670-690.
}
\seealso{
\code{\link{seBeta}}
}
\author{
Jeff Jones and Niels Waller
}
\keyword{Statistics}
