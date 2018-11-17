\name{RnpdMAP}
\alias{RnpdMAP}
\title{Generate Random NPD R matrices from a user-supplied population R
}
\description{
Generate a list of Random NPD (pseudo) R matrices with a user-defined fixed minimum eigenvalue from a user-supplied population R using the method of alternating projections.
}
\usage{
RnpdMAP(Rpop, 
        Lp = NULL, 
        NNegEigs = 1,
        NSmoothPosEigs = 4,
        NSubjects = NULL, 
        NSamples = 0, 
        MaxIts = 15000, 
        PRINT=FALSE, 
        Seed = NULL) 
}
\arguments{
  \item{Rpop}{input (PD or PSD) p x p Population correlation matrix.}
  \item{Lp}{desired minimum eigenvalue in the NPD matrices.}
  \item{NNegEigs}{number of eigenvalues < 0 in Rnpd.}
 \item{NSmoothPosEigs}{number of eigenvalues > 0 to smooth:  the smallest NSmoothPosEigs > 0be smoothed toward 0.}
 \item{NSubjects}{sample size (required when NSamples > 0) parameter used to generate sample correlation matrices. Default = NULL.}
 \item{NSamples}{generate NSamples sample R matrices. If NSamples = 0 the program will attempt to find Rnpd such that  ||Rpop - Rnpd||_2 is minimized.} 
 \item{MaxIts}{maximum number of projection iterations.}
 \item{PRINT}{(logical) If TURE the program will print the iteration history for Lp. Default = NULL.} 
 \item{Seed}{Optional seed for random number generation.}
 }
\value{
 \item{Rpop}{population (PD) correlation matrix.}
 \item{R}{sample correlation matrix.}
 \item{Rnpd}{NPD improper (pseudo) correlation matrix.}
 \item{Lp}{desired value of minimum eigenvalue.}
 \item{minEig}{observed value of minimum eigenvalue of Rnpd.}
\item{convergence}{0 = converged; 1 = not converged in MaxIts iterations of the alternating projections algorithm.}
\item{feasible}{logical) TRUE if max(abs(r_ij)) <= 1. If FALSE then one or more values in Rnpd > 1 in absolute value.}  
\item{Seed}{saved seed for random number generator.}
\item{prbs1}{vector probabilities used to generate eigenvalues < 0.}
\item{prbs2}{vector of probabilities used to smooth the smallest NSmoothPosEigs towards zero.}
}
\author{
Niels G. Waller
}
\examples{
library(MASS)

Nvar = 20
Nfac = 4
NSubj = 600
Seed = 123    

set.seed(Seed)

## Generate a vector of classical item difficulties
p <- runif(Nvar)

cat("\nClassical Item Difficulties:\n")

print(rbind(1:Nvar,round(p,2)) )

summary(p)


## Convert item difficulties to quantiles
b <- qnorm(p)


## fnc to compute root mean squared standard deviation
RMSD <- function(A, B){
  sqrt(mean( ( A[lower.tri(A, diag = FALSE)] - B[lower.tri(B, diag = FALSE)] )^2))
}


## Generate vector of eigenvalues with clear factor structure
  L <- eigGen(nDimensions = Nvar, 
            nMajorFactors = Nfac, 
            PrcntMajor = .60, 
            threshold  = .50)
          

## Generate a population R matrix with the eigenvalues in L
  Rpop <- rGivens(eigs = L)$R
  
## Generate continuous data that will reproduce Rpop (exactly)
  X <- mvrnorm(n = NSubj, mu = rep(0, Nvar), 
               Sigma = Rpop, empirical = TRUE)
               
if( any(colSums(X) == 0) ){
  stop("One or more variables have zero variance. Generate a new data set.")               
 }
 
## Cut X at thresholds given in b to produce binary data U
  U <- matrix(0, nrow(X), ncol(X))
  for(j in 1:Nvar){
    U[X[,j] <= b[j],j] <- 1
  }
  
## Compute tetrachoric correlations
  Rtet <- tetcor(U, Smooth = FALSE, PRINT = TRUE)$r
  # Calculate eigenvalues of tetrachoric R matrix
  Ltet <- eigen(Rtet)$values
  
  if(Ltet[Nvar] >= 0) stop("Rtet is P(S)D")
  
## Simulate NPD R matrix with minimum eigenvalue equal to 
  # min(Ltet)
  out <- RnpdMAP(Rpop, 
               Lp = Ltet[Nvar], 
               NNegEigs = Nvar/5,
               NSmoothPosEigs = Nvar/5, 
               NSubjects = 150, 
               NSamples = 1, 
               MaxIts = 15000, 
               PRINT = FALSE, 
               Seed = Seed) 

## RLp is a NPD pseudo R matrix with min eigenvalue = min(Ltet)
  RLp <- out[[1]]$Rnpd

## Calculate eigenvalues of simulated NPD R matrix (Rnpd)
  Lnpd <- eigen(RLp, only.values = TRUE)$values
  
## Scree plots for observed and simulated NPD R matrices.  
  ytop <- max(c(L,Lnpd,Ltet))
  pointSize = .8
  plot(1:Nvar, L, typ = "b", col = "darkgrey", lwd=3, 
       lty=1, 
       main = 
       "Eigenvalues of Rpop, Tet R, and Sim Tet R:
       \nSimulated vs Observed npd Tetrachoric R Matrices",
       ylim = c(-1, ytop),
       xlab = "Dimensions", 
       ylab = "Eigenvalues",
       cex = pointSize,cex.main = 1.2)
  points(1:Nvar, Lnpd, typ="b", 
         col = "red", lwd = 3, lty=2, cex=pointSize)
  points(1:Nvar, Ltet, typ="b", 
         col = "darkgreen", lwd = 3, lty = 3, cex= pointSize)
 
  legend("topright", 
         legend = c("eigs Rpop", "eigs Sim Rnpd", "eigs Emp Rnpd"), 
         col = c("darkgrey", "red","darkgreen"), 
         lty = c(1,2,3), 
         lwd = c(4,4,4), cex = 1.5)
  
  abline(h = 0, col = "grey", lty = 2, lwd = 4)
 
  cat("\nRMSD(Rpop, Rtet) = ", round(rmsd(Rpop, Rtet), 3))
  cat("\nRMSD(Rpop, RLp) = ",  round(rmsd(Rpop, RLp),  3))
}

\keyword{fungible}

