% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/faLocalMin.R
\name{faLocalMin}
\alias{faLocalMin}
\title{Investigate local minima in  faMain objects}
\usage{
faLocalMin(fout, Set = 1, HPthreshold = 0.1, digits = 5, PrintLevel = 1)
}
\arguments{
\item{fout}{(Object from class  \code{faMain}).}

\item{Set}{(Integer) The index of the solution set (i.e., the collection of 
rotated factor patterns with a common complexity value) from an 
\code{faMain} object.}

\item{HPthreshold}{(Scalar) A number between [0, 1] that defines the 
hyperplane threshold. Factor pattern elements below \code{HPthreshold} in absolute 
value are counted in the hyperplane count.}

\item{digits}{(Integer) Specifies the  number of significant 
digits in the printed output. Default \code{digits = 5}.}

\item{PrintLevel}{(Integer) Determines the level of printed output.
PrintLevel = 
\itemize{
 \item \strong{0}: No output is printed. 
 \item \strong{1}: Print output for the six most discrepant pairs of 
 rotated factor patterns.
 \item \strong{2}: Print output for all  pairs of rotated factor patterns.
}}
}
\value{
\code{faLocalMin} function will produce the following output.
\itemize{
  \item \strong{rmsdTable}: (Matrix) A table of \code{RMSD} values for each  pair of 
      rotated factor patterns in  solution set \code{Set}.
  \item \strong{Set}: (Integer) The index of the user-specified solution set.
  \item \strong{complexity.val} (Numeric): The common complexity value for all members 
  in the user-specified solution set.
  \item \strong{HPcount}: (Integer) The hyperplane count for each factor pattern in the solution set.
  }
}
\description{
Compute pairwise root mean squared deviations (RMSD) 
among rotated factor patterns in an \code{faMain} object. 
Prior to computing the RMSD values, each pair of solutions is aligned to 
the first member of the pair.  Alignment is accomplished using the 
Hungarian algorithm as described in \code{faAlign}.
}
\details{
Compute pairwise RMSD values among rotated factor patterns from 
an \code{faMain} object.
}
\examples{
\dontrun{
  ## Generate Population Model and Monte Carlo Samples ####
  sout <- simFA(Model = list(NFac = 5,
                          NItemPerFac = 5,
                           Model = "orthogonal"),
              Loadings = list(FacLoadDist = "fixed",
                              FacLoadRange = .8),
              MonteCarlo = list(NSamples = 100, 
                                SampleSize = 500),
              Seed = 655342)

  ## Population EFA loadings
  (True_A <- sout$loadings)

  ## Population Phi matrix
  sout$Phi

  ## Compute EFA on Sample 67 ####
  fout <- faMain (R = sout$Monte$MCData[[67]],
                numFactors = 5,
                targetMatrix = sout$loadings,
                facMethod = "fals",
                rotate= "cfT",
                rotateControl = list(numberStarts = 50,
                                     standardize="CM",
                                     kappa = 1/25),
                Seed=3366805)

  ## Summarize output from faMain
  summary(fout, Set = 1, DiagnosticsLevel = 2, digits=4)

  ## Investigate Local Solutions
  LMout <- faLocalMin(fout, 
                    Set = 1,
                    HPthreshold = .15,
                    digits= 5, 
                    PrintLevel = 1)
                    
  ## Print hyperplane count for each factor pattern 
  ## in the solution set
  LMout$HPcount
  }
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
Niels Waller
}
\concept{Factor Analysis Routines}
\keyword{Statistics}
\keyword{fungible}
