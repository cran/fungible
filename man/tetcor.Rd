\name{tetcor}
\alias{tetcor}
\title{Compute ML Tetrachoric Correlations}
\description{Compute ML tetrachoric correlations with optional bias correction and smoothing.}
\usage{
tetcor(X, y = NULL, BiasCorrect, stderror, Smooth = TRUE, max.iter)
}
\arguments{
  \item{X}{Either a matrix or vector of (0/1) binary data.}
  \item{y}{An optional(if X is a matrix) vector of (0/1) binary data.}
  \item{BiasCorrect}{A logical that determines whether bias correction (Brown & Benedetti, 1977) is performed. Default = TRUE.}
  \item{stderror}{A logical that determines whether standard errors are calulated.  Default = FALSE.}
  \item{Smooth}{A logical which determines whether the tetrachoric correlation matrix should be smoothed. 
  	            A smoothed matrix is always positive definite.}
  \item{max.iter}{Maximum number of iterations. Default = 50.}
}
\value{
 If stderror = FALSE, \code{tetcor} returns a matrix of tetrachoric correlations. If \code{stderror} = TRUE then \code{tetcor} returns a list the first component of which is a matrix of tetrachoric correlations and the second component is a matrix of standard errors (see Hamdan, 1970). 
}
\author{Niels Waller}
\references{
Brown, M. B. & Benedetti, J. K. (1977). On the mean and variance of the tetrachoric correlation coefficient. \emph{Psychometrika, 42}, 347--355.                              

Divgi, D. R. (1979) Calculation of the tetrachoric correlation coefficient. \emph{Psychometrika, 44}, 169-172.

Hamdan, M. A. (1970). The equivalence of tetrachoric and maximum likelihood estimates of rho in 2 by 2 tables. \emph{Biometrika, 57}, 212-215.                         
}

\examples{
## generate bivariate normal data
library(MASS)
set.seed(123)
rho <- .85
xy <- mvrnorm(100000, mu = c(0,0), Sigma = matrix(c(1, rho, rho, 1), ncol = 2))

# dichotomize at difficulty values
p1 <- .7
p2 <- .1
xy[,1] <- xy[,1] < qnorm(p1) 
xy[,2] <- xy[,2] < qnorm(p2)

print( apply(xy,2,mean), digits = 2)
#[1] 0.700 0.099

tetcor(X = xy, BiasCorrect = TRUE, 
       stderror = TRUE, Smooth = TRUE, max.iter = 5000)

# $r
# [,1]      [,2]
# [1,] 1.0000000 0.8552535
# [2,] 0.8552535 1.0000000
# 
# $se
# [,1]       [,2]
# [1,] 1.00000000 0.01458171
# [2,] 0.01458171 1.00000000
# 
# $Warnings
# list()


}
\keyword{Statistics}
