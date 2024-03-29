% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/faMain.R
\name{faMain}
\alias{faMain}
\title{Automatic Factor Rotation from Random Configurations with Bootstrap Standard Errors}
\usage{
faMain(
  X = NULL,
  R = NULL,
  n = NULL,
  numFactors = NULL,
  facMethod = "fals",
  urLoadings = NULL,
  rotate = "oblimin",
  targetMatrix = NULL,
  bootstrapSE = FALSE,
  numBoot = 1000,
  CILevel = 0.95,
  Seed = 1,
  digits = NULL,
  faControl = NULL,
  rotateControl = NULL,
  ...
)
}
\arguments{
\item{X}{(Matrix) A raw data matrix (or data frame).}

\item{R}{(Matrix) A correlation matrix.}

\item{n}{(Numeric) Sample size associated with the correlation matrix. Defaults to n = NULL.}

\item{numFactors}{(Numeric) The number of factors to extract for subsequent rotation.}

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

\item{urLoadings}{(Matrix) An unrotated factor-structure matrix to be rotated.}

\item{rotate}{(Character) Designate which rotation algorithm to apply. The 
following are available rotation options: "oblimin", "quartimin", "targetT", 
"targetQ", "oblimax", "entropy", "quartimax", "varimax", "simplimax", 
"bentlerT", "bentlerQ", "tandemI", "tandemII", "geominT", "geominQ", "cfT", 
"cfQ", "infomaxT", "infomaxQ", "mccammon", "bifactorT", "bifactorQ", and 
"none". Defaults to rotate = "oblimin". See \pkg{GPArotation} package for more 
details. Note that rotations ending in "T" and "Q" represent orthogonal and 
oblique rotations, respectively.}

\item{targetMatrix}{(Matrix) This argument serves two functions. First, if a 
user has requested either a "targetT" or "targetQ' rotation, then 
 the target matrix is used to conduct a fully or partially
specified target rotation. In the latter case,  freely estimated factor 
loadings are designated by "NA" values and rotation will be conducted using  
Browne's (1972a, 1972b, 2001) method for a partially-specified 
target rotation. Second, if any other rotation option is chosen then all 
rotated loadings matrices (and assorted output) will be aligned 
(but not rotated) with the target solution.}

\item{bootstrapSE}{(Logical) Computes bootstrap standard errors. All bootstrap 
samples are aligned to the global minimum solution. Defaults to 
bootstrapSE = FALSE (no standard errors).}

\item{numBoot}{(Numeric) The number bootstraps. Defaults to numBoot = 1000.}

\item{CILevel}{(Numeric) The confidence level (between 0 and 1) of the bootstrap 
confidence interval. Defaults to CILevel = .95.}

\item{Seed}{(Numeric) Starting seed for reproducible bootstrap results and factor rotations. 
Defaults to Seed = 1.}

\item{digits}{(Numeric) Rounds the values to the specified number of decimal 
places. Defaults to digits = NULL (no rounding).}

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

\item{...}{Values to be passed to the \code{\link[stats]{cor}} function.
\itemize{
  \item \strong{use}: (Character) A character string giving a method for 
  computing correlations in the presence of missing values: "everything" 
  (the default), "all.obs", "complete.obs", "na.or.complete", or 
  "pairwise.complete.obs".
  \item \strong{method}: (Character) A character string indicating which 
  correlation coefficient is to be computed: "pearson" (the default), 
  "kendall", or "spearman". 
  \item \strong{na.rm}: (Logical) Should missing values be removed (TRUE) 
  or not (FALSE)?
}}
}
\value{
The \code{faMain} function will produce a lot of output in addition 
to the rotated factor pattern matrix and the factor correlations.
\itemize{
  \item \strong{R}: (Matrix) Returns the correlation matrix, useful when raw data are supplied.
  \item \strong{loadings}: (Matrix) The rotated factor solution with the 
  lowest evaluated discrepancy function. This solution has the lowest 
  discrepancy function \emph{of the examined random starting configurations}. 
  It is not guaranteed to find the "true" global minimum. Note that multiple
   (or even all) local solutions can have the same discrepancy functions.
  \item \strong{Phi}: (Matrix) The factor correlations of the rotated factor 
  solution with the lowest evaluated discrepancy function (see Details).
  \item \strong{facIndeterminacy}: (Vector) A vector (with length equal to the number of factors)
  containing Guttman's (1955) index of factor indeterminacy for each factor. 
  \item \strong{h2}: (Vector) The vector of final communality estimates. 
  \item \strong{loadingsSE}: (Matrix) The matrix of factor-loading standard 
  errors across the bootstrapped factor solutions. Each matrix element is 
  the standard deviation of all bootstrapped factor loadings for that element position.
  \item \strong{CILevel} (Numeric) The user-defined confidence level (between 0 and 1) of the bootstrap 
   confidence interval. Defaults to CILevel = .95.
  \item \strong{loadingsCIupper}: (Matrix) Contains the upper confidence 
  interval of the bootstrapped factor loadings matrix. The confidence 
  interval width is specified by the user.
  \item \strong{loadingsCIlower}: (Matrix) Contains the lower confidence 
  interval of the bootstrapped factor loadings matrix. The confidence 
  interval width is specified by the user.
  \item \strong{PhiSE}: (Matrix) The matrix of factor correlation standard 
  errors across the bootstrapped factor solutions. Each matrix element is 
  the standard deviation of all bootstrapped factor correlations for that element position.
  \item \strong{PhiCIupper}: (Matrix) Contains the upper confidence interval 
  of the bootstrapped factor correlation matrix. The confidence interval 
  width is specified by the user.
  \item \strong{PhiCIlower}: (Matrix) Contains the lower confidence interval 
  of the bootstrapped factor correlation matrix. The confidence interval 
  width is specified by the user.
  \item \strong{facIndeterminacySE}: (Matrix) A row vector containing the 
  standard errors of Guttman's (1955) factor indeterminacy indices across the 
  bootstrap factor solutions. 
  \item \strong{localSolutions}: (List) A list containing all local solutions 
  in ascending order of their factor loadings, rotation complexity values (i.e., the first solution 
  is the "global" minimum). Each solution returns the 
  \itemize{
     \item \strong{loadings}: (Matrix) the factor loadings, 
     \item \strong{Phi}: (Matrix) factor correlations, 
     \item \strong{RotationComplexityValue}: (Numeric) the complexity value of the rotation algorithm, 
     \item \strong{facIndeterminacy}: (Vector) A vector of factor indeterminacy indices for each common factor, and 
     \item \strong{RotationConverged}: (Logical) convergence status of the rotation algorithm. 
     }
  \item \strong{numLocalSets} (Numeric) How many sets of local solutions
   with the same discrepancy value were obtained. 
  \item \strong{localSolutionSets}: (List) A list containing the sets of 
  unique local minima solutions. There is one list element for every unique 
  local solution that includes (a) the factor loadings matrix, (b) the factor 
  correlation matrix (if estimated), and (c) the discrepancy value of the rotation algorithm. 
  \item \strong{loadingsArray}: (Array) Contains an array of all bootstrapped 
  factor loadings. The dimensions are factor indicators, factors, and the 
  number of bootstrapped samples (representing the row, column, and depth, respectively).
  \item \strong{PhiArray}: (Array) Contains an array of all bootstrapped 
  factor correlations. The dimension are the number of factors, the number 
  of factors, and the number of bootstrapped samples (representing the row,
   column, and depth, respectively).
  \item \strong{facIndeterminacyArray}: (Array) Contains an array of all 
  bootstrap factor indeterminacy indices. The dimensions are 1, the number 
  of factors, and the number of bootstrap samples (representing the row, 
  column, and depth order, respectively).
  \item \strong{faControl}: (List) A list of the control parameters passed 
  to the factor extraction (\code{\link{faX}}) function.
  \item \strong{faFit}: (List) A list of additional output from the factor
  extraction routines. 
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
    \item \strong{convergedX}: (Logical) TRUE if the factor 
    \strong{extraction} routine converged. 
    \item \strong{convergedR}: (Logical) TRUE if the factor \strong{rotation} 
    routine converged (for the local solution with the minimum discrepancy 
    value).
  }
  \item \strong{rotateControl}: (List) A list of the control parameters 
  passed to the rotation algorithm.
  \item \strong{unSpunSolution}: (List) A list of output parameters (e.g., loadings, Phi, etc) from 
  the rotated solution that was obtained by rotating directly from the unrotated (i.e., unspun) common factor orientation. 
  \item \strong{targetMatrix} (Matrix) The input target matrix if supplied by the user.
  \item \strong{Call}: (call) A copy of the function call.
}
}
\description{
This function conducts factor rotations (using the \pkg{GPArotation} package) 
from a user-specified number of random (orthogonal) starting configurations. 
Based on the resulting complexity function value, the function determines the 
number of local minima and, among these local solutions, will find the 
"global minimum" (i.e., the minimized complexity value from the finite 
number of solutions). See Details below for an elaboration on the global 
minimum. This function can also return bootstrap standard errors of the factor solution.
}
\details{
\itemize{
  \item \strong{Global Minimum}: This function uses several random starting 
  configurations for factor rotations in an attempt to find the global 
  minimum solution. However, this function is not guaranteed to find the 
  global minimum. Furthermore, the global minimum solution need not be 
  more psychologically interpretable than any of the local solutions (cf. 
  Rozeboom, 1992). As is recommended, our function returns all local 
  solutions so users can make their own judgements.
  \item \strong{Finding clusters of local minima}: We find local-solution sets by sorting the rounded 
  rotation complexity values (to the number of  digits specified in the \code{epsilon} 
  argument of the \code{rotateControl} list) into sets with equivalent values. For example, 
  by default \code{epsilon = 1e-5.} will only evaluate the complexity 
  values to five significant digits. Any differences beyond that value will not effect the final sorting. 
}
}
\examples{
## Example 1

## Generate an oblique factor model
lambda <- matrix(c(.41, .00, .00,
                   .45, .00, .00,
                   .53, .00, .00,
                   .00, .66, .00,
                   .00, .38, .00,
                   .00, .66, .00,
                   .00, .00, .68,
                   .00, .00, .56,
                   .00, .00, .55),
                 nrow = 9, ncol = 3, byrow = TRUE)

## Generate factor correlation matrix
Phi <- matrix(.50, nrow = 3, ncol = 3)
diag(Phi) <- 1

## Model-implied correlation matrix
R <- lambda \%*\% Phi \%*\% t(lambda)
diag(R) <- 1

## Load the MASS package to create multivariate normal data
library(MASS)

## Generate raw data to perfectly reproduce R
X <- mvrnorm(Sigma = R, mu = rep(0, nrow(R)), empirical = TRUE, n = 300)

\dontrun{
## Execute 50 promax rotations from a least squares factor extraction
## Compute 100 bootstrap samples to compute standard errors and 
## 80 percent confidence intervals
Out1 <- faMain(X             = X,
               numFactors    = 3,
               facMethod     = "fals",
               rotate        = "promaxQ",
               bootstrapSE   = TRUE,
               numBoot       = 100,
               CILevel       = .80,
               faControl     = list(treatHeywood = TRUE),
               rotateControl = list(numberStarts = 2,  
                                    power        = 4,
                                    standardize  = "Kaiser"),
               digits        = 2)
Out1[c("loadings", "Phi")] 
}

## Example 2

## Load Thurstone's (in)famous box data
data(Thurstone, package = "GPArotation")

## Execute 5 oblimin rotations with Cureton-Mulaik standardization 
Out2 <- faMain(urLoadings    = box26,
               rotate        = "oblimin",
               bootstrapSE   = FALSE,
               rotateControl = list(numberStarts = 5,
                                    standardize  = "CM",
                                    gamma        = 0,
                                    epsilon      = 1e-6),
               digits        = 2)
               
Out2[c("loadings", "Phi")]     

## Example 3

## Factor matrix from Browne 1972
lambda <- matrix(c(.664,  .322, -.075,
                   .688,  .248,  .192,
                   .492,  .304,  .224,
                   .837, -.291,  .037,
                   .705, -.314,  .155,
                   .820, -.377, -.104,
                   .661,  .397,  .077,
                   .457,  .294, -.488,
                   .765,  .428,  .009), 
                 nrow = 9, ncol = 3, byrow = TRUE)   
                 
## Create partially-specified target matrix
Targ <- matrix(c(NA, 0,  NA,
                 NA, 0,  0,
                 NA, 0,  0,
                 NA, NA, NA,
                 NA, NA, 0,
                 NA, NA, NA,
                 .7, NA, NA,
                 0,  NA, NA,
                 .7, NA, NA), 
               nrow = 9, ncol = 3, byrow = TRUE)  
               
## Perform target rotation              
Out3 <- faMain(urLoadings   = lambda,
               rotate       = "targetT",
               targetMatrix = Targ,
               digits       = 3)$loadings
Out3
}
\references{
Browne, M. W.  (1972).  Oblique rotation to a partially specified 
target.  \emph{British Journal of Mathematical and Statistical Psychology, 25},(1), 
207-212.

Browne, M. W. (1972b). Orthogonal rotation to a partially specifed
target. \emph{British Journal of Statistical Psychology, 25},(1), 115-120.

Browne, M. W. (2001). An overview of analytic rotation in 
exploratory factor analysis. \emph{Multivariate Behavioral Research, 36}(1), 111-150.

Cureton, E. E., & Mulaik, S. A. (1975). The weighted varimax 
rotation and the promax rotation. \emph{Psychometrika, 40}(2), 183-195.

Guttman, L. (1955). The determinacy of factor score matrices with 
implications for five other basic problems of common factor theory. 
\emph{British Journal of Statistical Psychology, 8}(2), 65-81.

Jung, S. & Takane, Y.  (2008).  Regularized common factor analysis.  
\emph{New Trends in Psychometrics}, 141-149.

Mansolf, M., & Reise, S. P. (2016). Exploratory bifactor 
analysis: The Schmid-Leiman orthogonalization and Jennrich-Bentler 
analytic rotations. \emph{Multivariate Behavioral Research, 51}(5), 698-717.

Rozeboom, W. W. (1992). The glory of suboptimal factor rotation: 
Why local minima in analytic optimization of simple structure are more 
blessing than curse. \emph{Multivariate Behavioral Research, 27}(4), 585-599.

Zhang, G. (2014). Estimating standard errors in exploratory factor 
analysis. \emph{Multivariate Behavioral Research, 49}(4), 339-353.
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
  \item Casey Giordano (Giord023@umn.edu)
  \item The authors thank Allie Cooperman and Hoang
   Nguyen for their help implementing the standard error estimation and the 
   Cureton-Mulaik standardization procedure.
}
}
\concept{Factor Analysis Routines}
