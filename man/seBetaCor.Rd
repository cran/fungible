\name{seBetaCor}
\alias{seBetaCor}
\title{Standard Errors and CIs for Standardized Regression Coefficients from Correlations}
\description{
 Computes Normal Theory and ADF Standard Errors and CIs for Standardized Regression Coefficients from Correlations
}
\usage{
seBetaCor(R,rxy, Nobs, alpha=.05,  digits=3, covmat = 'normal')
}
\arguments{
 \item{R}{A p x p predictor correlation matrix.}
 \item{rxy}{A p x 1 vector of predictor-criterion correlations }
 \item{Nobs}{Number of observations.}
 \item{alpha}{Desired Type I error rate; default = .05.}
 \item{digits}{Number of significant digits to print; default = 3. }
 \item{covmat}{String = 'normal' (the default) or a (p+1)p/2 x (p+1)p/2 covariance matrix of correlations.  The default option computes an asymptotic covariance matrix under the assumption of multivariate normal data. Users can supply a covariance matrix under asymptotic distribution free (ADF) or elliptical distributions when available. }
}
\value{
\item{cov.Beta}{Covariance matrix of standardized regression coefficients.}
\item{se.Beta}{Vector of standard errors for the standardized regression coefficients.} 
\item{alpha}{Type-I error rate. }                           
\item{CI.Beta}{(1-alpha)\% confidence intervals for standardized 
     regression coefficients. }
}
\references{
Jones, J. A, and Waller, N. G. (2013). The Normal-Theory and asymptotic 
distribution-free (ADF) covariance matrix of standardized regression coefficients: Theoretical extensions and finite sample behavior.Technical Report (052913)[TR052913]

Nel, D.A.G. (1985). A matrix derivation of the asymptotic covariance matrix of sample correlation coefficients. \emph{Linear Algebra and its Applications, 67}, 137-145.

Yuan, K. and Chan, W. (2011). Biases and standard errors of standardized regression coefficients. \emph{Psychometrika}, 76(4), 670--690.
}
\author{Jeff Jones and Niels Waller}
\examples{
R <- matrix(c(1.0000, 0.3511, 0.3661,
	          0.3511, 1.0000, 0.4359,
	          0.3661, 0.4359, 1.0000), 3, 3) 

rxy <- c(0.5820, 0.6997, 0.7621)
Nobs <- 46
out <- seBetaCor(R = R, rxy = rxy, Nobs = Nobs) 

# 95\% CIs for Standardized Regression Coefficients: 
#
#        lbound estimate ubound
# beta_1  0.107    0.263  0.419
# beta_2  0.231    0.391  0.552
# beta_3  0.337    0.495  0.653
}

\keyword{Statistics}% at least one, from doc/KEYWORDS
