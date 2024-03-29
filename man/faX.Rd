% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/faX.R
\name{faX}
\alias{faX}
\title{Factor Extraction (faX) Routines}
\usage{
faX(R, n = NULL, numFactors = NULL, facMethod = "fals", faControl = NULL)
}
\arguments{
\item{R}{(Matrix) A correlation matrix used for factor extraction.}

\item{n}{(Numeric) Sample size associated with the correlation matrix. 
Defaults to n = NULL.}

\item{numFactors}{(Numeric) The number of factors to extract for subsequent rotation.}

\item{facMethod}{(Character) The method used for factor extraction. The 
supported options are "fals" for unweighted least squares, "faml" for maximum 
likelihood, "fapa" for iterated principal axis factoring, and "pca" for 
principal components analysis. The default method is "fals".
\itemize{
  \item \strong{"fals"}: Factors are extracted using the unweighted least 
  squares estimation procedure using the \code{\link{fals}} function.
  \item \strong{"faml"}: Factors are extracted using the maximum likelihood 
  estimation procedure using the \code{\link[stats]{factanal}} function.
  \item \strong{"faregLS"}: Factors are extracted using regularized 
  least squares factor analysis using the \code{\link{fareg}} function. 
  \item \strong{"faregML"}: Factors are extracted using regularized 
  maximum likelihood factor using the \code{\link{fareg}} function. 
  \item \strong{"fapa"}: Factors are extracted using the iterated principal 
  axis factoring estimation procedure using the \code{\link{fapa}} function.
  \item \strong{"pca"}: Principal components are extracted. 
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
This function returns a list of output relating to the extracted factor loadings.
\itemize{
  \item \strong{loadings}: (Matrix) An unrotated factor structure matrix.
  \item \strong{h2}: (Vector) Vector of final communality estimates.
  \item \strong{faFit}: (List) A list of additional factor extraction output.
  \itemize{
    \item \strong{facMethod}: (Character) The factor extraction routine.
    \item \strong{df}: (Numeric) Degrees of Freedom from the maximum 
    likelihood factor extraction routine.
    \item \strong{n}: (Numeric) Sample size associated with the correlation matrix.
    \item \strong{objectiveFunc}: (Numeric) The evaluated objective function for the 
    maximum likelihood factor extraction routine. 
    \item \strong{RMSEA}: (Numeric) Root mean squared error of approximation 
    from Steiger & Lind (1980). Note that bias correction is computed if the 
    sample size is provided.
    \item \strong{testStat}: (Numeric) The significance test statistic for the maximum 
    likelihood procedure. Cannot be computed unless a sample size is provided. 
    \item \strong{pValue}: (Numeric) The p value associated with the significance test 
    statistic for the maximum likelihood procedure. Cannot be computed unless 
    a sample size is provided. 
    \item \strong{gradient}: (Matrix) The solution gradient for the least squares factor 
    extraction routine. 
    \item \strong{maxAbsGradient}: (Numeric) The maximum absolute value of the 
    gradient at the least squares solution. 
    \item \strong{Heywood}: (Logical) TRUE if a Heywood case was produced.
    \item \strong{converged}: (Logical) TRUE if the least squares or 
    principal axis factor extraction routine converged. 
  }
}
}
\description{
This function can be used to extract an unrotated factor structure matrix 
using the following algorithms: (a) unweighted least squares ("fals"); 
(b) maximum likelihood ("faml"); (c) iterated principal axis factoring ("fapa");
and (d) principal components analysis ("pca").
}
\details{
\itemize{
  \item \strong{Initial communality estimate}: According to Widaman and 
  Herringer (1985), the initial communality estimate does not have much 
  bearing on the resulting solution \emph{when the a stringent convergence 
  criterion is used}. In their analyses, a convergence criterion of .001 
  (i.e., slightly less stringent than the default of 1e-4) is sufficiently 
  stringent to produce virtually identical communality estimates irrespective 
  of the initial estimate used. It should be noted that all four methods for 
  estimating the initial communality in Widaman and Herringer (1985) are the 
  exact same used in this function. Based on their findings, it is not 
  recommended to use a convergence criterion lower than 1e-3.
}
}
\examples{
## Generate an example factor structure matrix
lambda <- matrix(c(.62, .00, .00,
                   .54, .00, .00,
                   .41, .00, .00,
                   .00, .31, .00,
                   .00, .58, .00,
                   .00, .62, .00,
                   .00, .00, .38,
                   .00, .00, .43,
                   .00, .00, .37),
                 nrow = 9, ncol = 3, byrow = TRUE)

## Find the model implied correlation matrix
R <- lambda \%*\% t(lambda)
diag(R) <- 1

## Extract (principal axis) factors using the factExtract function
Out1 <- faX(R          = R,
            numFactors = 3,
            facMethod  = "fapa",
            faControl  = list(communality = "maxr",
                              epsilon     = 1e-4))

## Extract (least squares) factors using the factExtract function
Out2 <- faX(R          = R,
            numFactors = 3,
            facMethod  = "fals",
            faControl  = list(treatHeywood = TRUE))

}
\references{
Jung, S. & Takane, Y.  (2008).  Regularized common factor analysis.  
New trends in psychometrics, 141-149.

Steiger, J. H., & Lind, J. (1980). Paper presented at the annual 
meeting of the Psychometric Society. \emph{Statistically-based tests for the 
number of common factors.}

Widaman, K. F., & Herringer, L. G. (1985). Iterative least squares 
estimates of communality: Initial estimate need not affect stabilized value. 
\emph{Psychometrika, 50}(4), 469-477.
}
\seealso{
Other Factor Analysis Routines: 
\code{\link{BiFAD}()},
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
  \item Casey Giordano (Giord023@umn.edu)
  \item Niels G. Waller (nwaller@umn.edu)
}
}
\concept{Factor Analysis Routines}
