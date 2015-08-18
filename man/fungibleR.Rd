\name{fungibleR}
\alias{fungibleR}
\title{Generate Fungible Correlation Matrices}
\description{
Generate fungible correlation matrices. For a given vector of standardized regression coefficients, Beta,  and a user-define R-squared value, Rsq, find predictor correlation matrices, R, such that Beta' R Beta = Rsq. The size of the smallest eigenvalue (Lp) of R can be defined.  
}
\usage{
fungibleR(R, Beta, Lp = .00, eps = 1e-16)
}
\arguments{
 \item{R}{A p x p predictor correlation matrix.}
 \item{Beta}{A p x 1 vector of standardized regression coefficients.}
 \item{Lp}{Controls the size of the smallest eigenvalue of RstarLp.} 
 \item{eps}{Convergence criterion.}
}
\value{
  \item{R}{Any input correlation matrix that satisfies Beta' R Beta = Rsq.}
  \item{Beta}{Input vector of std reg coefficients.}
  \item{Rstar}{A random fungible correlation matrix.}
  \item{RstarLp}{A fungible correlation matrix with a fixed minimum eigenvalue
                (RstarLp can be PD, PSD, or ID).}
  \item{s}{Scaling constant for Rstar.}
  \item{sLp}{Scaling constant for RstarLp.}
  \item{Delta}{Vector in the null space of vecp(Beta Beta').}
  \item{Q}{Left null space of Beta.}
  \item{FrobNorm}{Frobenius norm ||R - Rstar||_F.}
  \item{FrobNormLp}{Frobenius norm ||R - RstarLp||_F given random Delta.}
  \item{converged}{An integer code. 0 indicates successful completion.}
}
\references{
Waller, N. (2015). Fungible Correlation Matrices: A New Tool for Evaluating Penalized Regression Models.
}
\author{Niels Waller}
\examples{
R <- matrix(.5, 5, 5)
diag(R) <- 1
Beta <- rep(.1, 5)

## change value of Lp to control the size of the smallest eigenvalue of RstarMax
out <- fungibleR(R, Beta, Lp = 0.12345678, eps = 1e-16)

## R
print(round(out$R,3)) 

## Rstar
print(round(out$Rstar,3)) 

## Eigen Rstar
print(round(eigen(out$Rstar)$values, 9)) 

## RstarLp
print(round(out$RstarLp, 3)) 

## Eigen RstarMax
print(eigen(out$RstarLp)$values, digits = 14) 

## Coefficient of determination (R^2)
t(Beta) \%*\% out$RstarLp \%*\% Beta 

if(out$converged) print("Falied to converge")
}
\keyword{fungible}% at least one, from doc/KEYWORDS
