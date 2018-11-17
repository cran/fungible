\name{smoothLG}
\alias{smoothLG}
\title{Smooth NPD to Nearest PSD or PD Matrix}
\description{Smoothing an indefinite matrix to a PSD matrix via theory described by Lurie and Goldberg}
\usage{ 
smoothLG(R, start.val = NULL, Wghts = NULL, PD = FALSE, 
	Penalty = 50000, eps=1e-07)
}
\arguments{
\item{R}{Indefinite Matrix.}
\item{start.val}{Optional vector of start values for Cholesky factor of S.}  
\item{Wghts}{An optional matrix of weights such that the objective function minimizes wij(rij - sij)^2, where wij is Wghts[i,j].}  
\item{PD}{Logical (default = FALSE). If PD = TRUE then the objective function  will smooth the least squares solution to insure Positive Definitness.}
\item{Penalty}{A scalar weight to scale the Lagrangian multiplier. Default = 50000.}
\item{eps}{A small value to add to zero eigenvalues if smoothed matrix must be PD. Default = 1e-07.}  
}
\value{
\item{RLG}{Lurie Goldberg smoothed matrix.}
\item{RKB}{Knol and Berger smoothed matrix.}
\item{convergence}{0 = converged solution, 1 = convergence failure.}
\item{start.val}{Vector of start.values.}
\item{gr}{Analytic gradient at solution.} 
\item{Penalty}{Scalar used to scale the Lagrange multiplier.}  
\item{PD}{User-supplied value of PD.}  
\item{Wghts}{Weights used to scale the squared euclidean distances.}
\item{eps}{Value added to zero eigenvalue to produce PD matrix.} 
}
\author{Niels Waller}
\examples{
data(BadRLG)

out<-smoothLG(R = BadRLG, Penalty = 50000)
cat("\nGradient at solution:", out$gr,"\n")
cat("\nNearest Correlation Matrix\n")
print( round(out$RLG,8) )

################################
##  Rousseeuw Molenbergh example
data(BadRRM)

out <- smoothLG(R = BadRRM, PD=TRUE)
cat("\nGradient at solution:", out$gr,"\n")
cat("\nNearest Correlation Matrix\n")
print( round(out$RLG,8) )

## Weights for the weighted solution
W <- matrix(c(1,  1, .5,
              1,  1,  1,
              .5,  1,  1), nrow = 3, ncol = 3)
tmp <- smoothLG(R = BadRRM,  PD = TRUE, eps=.001)
cat("\nGradient at solution:", out$gr,"\n")
cat("\nNearest Correlation Matrix\n")
print( round(out$RLG,8) )
print( eigen(out$RLG)$val )

## Rousseeuw Molenbergh 
## non symmetric matrix
T <- matrix(c(.8, -.9, -.9, 
            -1.2,  1.1, .3, 
             -.8, .4, .9),  nrow = 3, ncol = 3,byrow=TRUE)
out <- smoothLG(R = T,  PD = FALSE, eps=.001)

cat("\nGradient at solution:", out$gr,"\n")
cat("\nNearest Correlation Matrix\n")
print( round(out$RLG,8) )    
}
\keyword{statistics}% at least one, from doc/KEYWORDS

