% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smoothBY.R
\name{smoothBY}
\alias{smoothBY}
\title{Smooth an NPD R matrix to PD using the Bentler Yuan 2011 method}
\usage{
smoothBY(R, const = 0.98, eps = 0.001)
}
\arguments{
\item{R}{Indefinite Matrix.}

\item{const}{const is a user-defined parameter that is defined as k in
Bentler and Yuan (2011). If 0 < const < 1, then const is treated as a fixed
value. If const = 1 then the program will attempt to find the highest value
of const such that R is positive (semi) definite.}

\item{eps}{If const = 1 then the program will iteratively reduce const by
eps until either (a) the program converges or (b) const < = 0.}
}
\value{
\item{RBY}{smoothed correlation matrix.} \item{constant}{The final
value of const.} \item{convergence}{(Logical) a value of TRUE indicates that
the function converged.} \item{outStatus}{Convergence state for Rcsdp::csdp
function. \cr \cr 0: \cr \cr Success. Problem solved to full accuracy \cr
\cr 1: \cr \cr Success. Problem is primal infeasible \cr \cr 2: \cr \cr
Success. Problem is dual infeasible \cr \cr 3: \cr \cr Partial Success.
Solution found but full accuracy was not achieved \cr \cr 4: \cr \cr
Failure. Maximum number of iterations reached \cr \cr 5: \cr \cr Failure.
Stuck at edge of primal feasibility \cr \cr 6: \cr \cr Failure. Stuch at
edge of dual infeasibility \cr \cr 7: \cr \cr Failure. Lack of progress \cr
\cr 8:\cr \cr Failure. X or Z (or Newton system O) is singular \cr \cr 9:
\cr \cr Failure. Detected NaN or Inf values} \cr \item{glb}{Greatest lower
bound reliability estimates.} \item{eps}{Default value (eps = 1E-03) or
user-supplied value of eps.}
}
\description{
Smooth a NPD correlation matrix to PD using the Bentler and Yuan method.
}
\examples{

data(BadRBY)

out<-smoothBY(R = BadRBY, const = .98)
cat("\nSmoothed Correlation Matrix\n")
print( round(out$RBY,8) )
cat("\nEigenvalues of smoothed matrix\n")
print( eigen(out$RBY)$val  )

}
\references{
Bentler, P. M. & Yuan, K. H.  (2011).  Positive definiteness via
off-diagonal scaling of a symmetric indefinite matrix.  \emph{Psychometrika,
76}(1), 119--123.

Debelak, R. & Tran, U. S. (2013). Principal component analysis of smoothed
tetrachoric correlation matrices as a measure of dimensionality.
\emph{Educational and Psychological Measurement, 73}(1), 63--77.
}
\author{
Code modified from that reported in Debelak, R. & Tran, U. S.
(2011).
}
\keyword{statistics}
