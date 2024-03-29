% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BiFAD.R
\name{BiFAD}
\alias{BiFAD}
\title{Bifactor Analysis via Direct Schmid-Leiman (DSL) Transformations}
\usage{
BiFAD(
  R,
  B = NULL,
  numFactors = NULL,
  facMethod = "fals",
  rotate = "oblimin",
  salient = 0.25,
  rotateControl = NULL,
  faControl = NULL
)
}
\arguments{
\item{R}{(Matrix) A correlation matrix.}

\item{B}{(Matrix) Bifactor target matrix. If B is NULL the program will create an empirically defined target matrix.}

\item{numFactors}{(Numeric) The number of group factors to estimate.}

\item{facMethod}{(Character) The method used for factor extraction 
(\code{\link{faX}}). The supported options are "fals" for unweighted least 
squares, "faml" for maximum likelihood, "fapa" for iterated principal axis 
factoring, "faregLS" for regularized least squares,
"faregML" for regularized maximum likelihood, and "pca" for principal components 
 analysis. The default method  is "fals". 
\itemize{
  \item \strong{"fals"}: Factors are extracted using the unweighted least 
  squares estimation procedure using the \code{\link{fals}} function.
  \item \strong{"faml"}: Factors are extracted using the maximum likelihood 
  estimation procedure using the \code{\link[stats]{factanal}} function.
  \item \strong{"fapa"}: Factors are extracted using the iterated principal 
  axis factoring estimation procedure using the \code{\link{fapa}} function.
  \item \strong{"faregLS"}: Factors are extracted using regularized 
  least squares factor analysis using the \code{\link{fareg}} function. 
  \item \strong{"faregML"}: Factors are extracted using regularized 
  maximum likelihood factor using the \code{\link{fareg}} function. 
  \item \strong{"pca"}: Principal components are extracted. 
}}

\item{rotate}{(Character) Designate which rotation algorithm to apply. See the \code{\link{faMain}} function for more details about possible rotations. An oblimin rotation is the default.}

\item{salient}{(Numeric) Threshold value for creating an empirical target matrix.}

\item{rotateControl}{(List) A list of control values to pass to the factor rotation algorithms.
\itemize{
  \item \strong{numberStarts}: (Numeric) The number of random (orthogonal) 
  starting configurations for the chosen rotation method (e.g., oblimin). The first
  rotation will always commence from the unrotated factors orientation.
  Defaults to numberStarts = 10. 
  \item \strong{gamma}: (Numeric) This is a tuning parameter (between 0 
  and 1, inclusive) for an oblimin rotation.  See the \pkg{GPArotation} 
  library's oblimin documentation for more details. Defaults to gamma = 0 
  (i.e., a quartimin rotation).
  \item \strong{delta}: (Numeric) This is a tuning parameter for the geomin
   rotation. It adds a small number (default = .01) to the squared factor 
   loadings before computing the geometric means in the discrepancy function.
  \item \strong{kappa}: (Numeric) The main parameterization of the 
  Crawford-Ferguson (CF) rotations (i.e., "cfT" and "cfQ" for orthogonal and 
  oblique CF rotation, respectively). Defaults to kappa = 0. 
  \item \strong{k}: (Numeric) A specific parameter of the simplimax rotation. 
  Defaults to k = the number of observed variables.
  \item \strong{standardize}: (Character) The standardization routine used 
  on the unrotated factor structure. The three options are "none", "Kaiser", 
  and "CM". Defaults to standardize = "none". 
  \itemize{
    \item \strong{"none"}: No standardization is applied to the unrotated 
    factor structure. 
    \item \strong{"Kaiser"}: Use a factor structure matrix that has been 
    normed by Kaiser's method (i.e., normalize all rows to have a unit length).
    \item \strong{"CM"}: Use a factor structure matrix that has been normed
     by the Cureton-Mulaik method.
  }
  \item \strong{epsilon}: (Numeric) The rotational convergence criterion to 
  use. Defaults to epsilon = 1e-5.
  \item \strong{power}: (Numeric) Raise factor loadings the the n-th power 
  in the \code{\link{promaxQ}} rotation. Defaults to power = 4.
  \item \strong{maxItr}: (Numeric) The maximum number of iterations for the 
  rotation algorithm. Defaults to maxItr = 15000.
}}

\item{faControl}{(List) A list of optional parameters passed to the factor 
extraction (\code{\link{faX}}) function.
\itemize{
  \item \strong{treatHeywood}: (Logical) In \code{fals}, if treatHeywood is 
  true, a penalized least squares function is used to bound the communality 
  estimates below 1.0. Defaults to treatHeywood = TRUE.
  \item \strong{nStart}: (Numeric) The number of starting values to be tried 
  in \code{faml}. Defaults to nStart = 10.
  \item \strong{start}: (Matrix) NULL or a matrix of starting values, each column 
  giving an initial set of uniquenesses. Defaults to start = NULL. 
  \item \strong{maxCommunality}: (Numeric) In \code{faml}, set the maximum 
  communality value for the estimated solution. Defaults to maxCommunality = .995.
  \item \strong{epsilon}: (Numeric) In \code{fapa}, the numeric threshold 
  designating when the algorithm has converged. Defaults to epsilon = 1e-4.
  \item \strong{communality}: (Character) The method used to estimate the 
  initial communality values in \code{fapa}. Defaults to communality = 'SMC'.
  \itemize{
    \item \strong{"SMC"}: Initial communalities are estimated by taking the 
    squared multiple correlations of each indicator after regressing the 
    indicator on the remaining variables.
    \item \strong{"maxr"}: Initial communalities equal the largest 
    (absolute value) correlation in each column of the correlation matrix.
    \item \strong{"unity"}: Initial communalities equal 1.0 for all variables.
  }
  \item \strong{maxItr}: (Numeric) In \code{fapa}, the maximum number of 
  iterations to reach convergence. Defaults to maxItr = 15,000.
}}
}
\value{
The following output are returned in addition to the estimated Direct Schmid-Leiman bifactor solution.

\itemize{
  \item \strong{B}: (Matrix) The target matrix used for the Procrustes rotation.
  \item \strong{BstarSL}: (Matrix) The resulting (rank-deficient) matrix of Direct Schmid-Leiman factor loadings.
  \item \strong{BstarFR}: (Matrix) The resulting (full-rank) matrix of Direct Bifactor factor loadings.
  \item \strong{rmsrSL}: (Scalar) The root mean squared residual (rmsr) between the known B matrix and the estimated (rank-deficient) Direct Schmid-Leiman rotation. If the B target matrix is empirically generated, this value is NULL.
  \item \strong{rmsrFR}: (Scalar) The root mean squared residual (rmsr) between the known B matrix and the estimated (full-rank) Direct Bifactor rotation. If the B target matrix is empirically generated, this value is NULL.
}
}
\description{
This function estimates the (rank-deficient) Direct Schmid-Leiman (DSL) bifactor solution as well as the (full-rank) Direct Bifactor (DBF) solution.
}
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

out.ex1 <- BiFAD(R          = Rex1,
                 B          = NULL,
                 numFactors = 4,
                 facMethod  = "fals",
                 rotate     = "oblimin",
                 salient    = .25)

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

out.ex2 <- BiFAD(R          = Rex1,
                 B          = Bpattern,
                 numFactors = NULL,
                 facMethod  = "fals",
                 rotate     = "oblimin",
                 salient    = .25)

cat("\nRank Deficient Bifactor Solution:\n")
print( round(out.ex2$BstarSL, 2) )

cat("\nFull Rank Bifactor Solution:\n")
print( round(out.ex2$BstarFR, 2) )

}
\references{
\itemize{
  \item Giordano, C. & Waller, N. G. (under review). Recovering bifactor models: A comparison of seven methods.
  \item Mansolf, M., & Reise, S. P. (2016). Exploratory bifactor analysis: The Schmid-Leiman orthogonalization and Jennrich-Bentler analytic rotations. \emph{Multivariate Behavioral Research, 51}(5), 698-717.
  \item Waller, N. G. (2018). Direct Schmid Leiman transformations and rank deficient loadings matrices. \emph{Psychometrika, 83}, 858-870.
}
}
\seealso{
Other Factor Analysis Routines: 
\code{\link{Box26}},
\code{\link{GenerateBoxData}()},
\code{\link{Ledermann}()},
\code{\link{SLi}()},
\code{\link{SchmidLeiman}()},
\code{\link{faAlign}()},
\code{\link{faEKC}()},
\code{\link{faIB}()},
\code{\link{faLocalMin}()},
\code{\link{faMB}()},
\code{\link{faMain}()},
\code{\link{faScores}()},
\code{\link{faSort}()},
\code{\link{faStandardize}()},
\code{\link{faX}()},
\code{\link{fals}()},
\code{\link{fapa}()},
\code{\link{fareg}()},
\code{\link{fsIndeterminacy}()},
\code{\link{orderFactors}()},
\code{\link{print.faMB}()},
\code{\link{print.faMain}()},
\code{\link{promaxQ}()},
\code{\link{summary.faMB}()},
\code{\link{summary.faMain}()}
}
\author{
\itemize{
  \item Niels G. Waller (nwaller@umn.edu)
}
}
\concept{Factor Analysis Routines}
