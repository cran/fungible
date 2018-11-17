\name{BiFAD}
\alias{BiFAD}
\title{Bifactor Analysis via Direct Schmid Leiman Transformations}
\description{
Bifactor Analysis via Direct Schmid Leiman Transformations 
}
\usage{
         BiFAD(R, B = NULL, nGroup = NULL, 
                  factorMethod = "minres", 
                  rotation="oblimin",  salient = .25, 
                  maxitFA = 5000, 
                  maxitRotate = 5000,
                  gamma = 0)
}

\arguments{
  \item{R}{Input correlation matrix.}
  \item{B}{Bifactor target matrix. If B=NULL the program will create an empirically defined target matrix.}
  \item{nGroup}{Number of group factors in bifactor solution.}
  \item{factorMethod}{Factor extraction method. Options include:  
  minres (minimum residual), ml (maximum likelihood), 
  pa (principal axis), gls (generalized least squares).}
 \item{rotation}{Factor rotation method. Current options include:
 oblimin, geominQ, quartimin, promax.}
 \item{salient}{Threshold value for creating an empirical target matrix.}
 \item{maxitFA}{Maximum iterations for the factor extraction method.} 
 \item{maxitRotate}{Maximum iterations for the gradient pursuit rotation algorithm.} 
 \item{gamma}{Optional tuning parameter for oblimin rotation.} 
}

\value{
  \item{B}{User defined or empirically generated target matrix.}
  \item{BstarSL}{Direct S-L solution.}
  \item{BstarFR}{Direct full rank bifactor solution.}
  \item{rmsrSL}{Root mean squared residual of (B - BstarSL).}
  \item{rmsrFR}{Root mean squared residual of (B - BstarFR).}
}
\references{
Waller, N. G. (2018). Direct Schmid Leiman transformations and rank deficient loadings matrices. \emph{Psychometrika}, 83, 858-870.
}
\author{Niels Waller}
\examples{
cat("\nExample 1:\nEmpirical Target Matrix:\n")
# Mansolf and Reise Table 2 Example
Btrue <- matrix(c(.48, .40,  0,   0,   0,
                  .51, .35,  0,   0,   0,
                  .67, .62,  0,   0,   0,
                  .34, .55,  0,   0,   0,
                  .44,  0, .45,   0,   0,
                  .40,  0, .48,   0,   0,
                  .32,  0, .70,   0,   0,
                  .45,  0, .54,   0,   0,
                  .55,  0,   0, .43,   0,
                  .33,  0,   0, .33,   0,
                  .52,  0,   0, .51,   0,
                  .35,  0,   0, .69,   0,
                  .32,  0,   0,   0, .65,
                  .66,  0,   0,   0, .51,
                  .68,  0,   0,   0, .39,
                  .32,  0,   0,   0, .56), 16, 5, byrow=TRUE)

Rex1 <- Btrue \%*\% t(Btrue)
diag(Rex1) <- 1

out.ex1 <- BiFAD(R = Rex1, B = NULL, nGroup = 4, 
                 factorMethod = "minres", 
                 rotation="oblimin",  salient = .25, 
                 maxitFA = 5000, 
                 maxitRotate = 5000,
                 gamma = 0)

cat("\nRank Deficient Bifactor Solution:\n")
print( round(out.ex1$BstarSL, 2) )

cat("\nFull Rank Bifactor Solution:\n")
print( round(out.ex1$BstarFR, 2) )

cat("\nExample 2:\nUser Defined Target Matrix:\n")

Bpattern <- matrix(c( 1,  1,  0,   0,   0,
                      1,  1,  0,   0,   0,
                      1,  1,  0,   0,   0,
                      1,  1,  0,   0,   0,
                      1,  0,  1,   0,   0,
                      1,  0,  1,   0,   0,
                      1,  0,  1,   0,   0,
                      1,  0,  1,   0,   0,
                      1,  0,   0,  1,   0,
                      1,  0,   0,  1,   0,
                      1,  0,   0,  1,   0,
                      1,  0,   0,  1,   0,
                      1,  0,   0,   0,  1,
                      1,  0,   0,   0,  1,
                      1,  0,   0,   0,  1,
                      1,  0,   0,   0,  1), 16, 5, byrow=TRUE)

out.ex2 <- BiFAD(R = Rex1, B = Bpattern, nGroup = NULL, 
                 factorMethod = "minres", 
                 rotation="oblimin",  salient = .25, 
                 maxitFA = 5000, 
                 maxitRotate = 5000,
                 gamma = 0)

cat("\nRank Deficient Bifactor Solution:\n")
print( round(out.ex2$BstarSL, 2) )

cat("\nFull Rank Bifactor Solution:\n")
print( round(out.ex2$BstarFR, 2) )
}
\keyword{BiFAD}% at least one, from doc/KEYWORDS

