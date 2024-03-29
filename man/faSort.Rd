% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/faSort.R
\name{faSort}
\alias{faSort}
\title{Sort a factor loadings matrix}
\usage{
faSort(fmat, phi = NULL, BiFactor = FALSE, salient = 0.25, reflect = TRUE)
}
\arguments{
\item{fmat}{factor loadings (pattern or structure) matrix.}

\item{phi}{factor correlation matrix. Default = NULL. If reflect = TRUE then
phi will be corrected to match the new factor orientations.}

\item{BiFactor}{(logical) Is the solution a bifactor model?}

\item{salient}{factor markers with loadings >= abs(salient) will be saved in
the markers list. Note that a variable can be a marker of more than one
factor.}

\item{reflect}{(logical) if reflect = TRUE then the factors will be
reflected such that salient loadings are mostly positive. Default \code{Reflect = TRUE}.}
}
\value{
\item{loadings}{sorted factor loadings matrix.} 
   \item{phi}{reflected factor correlation matrix when phi is given as an argument.}
   \item{markers}{A list of factor specific markers with loadings >=
abs(salient). Markers are sorted by the absolute value of the salient factor
loadings.} \item{sortOrder}{sorted row numbers.} 
   \item{SEmat}{The SEmat is a
so-called Start-End matrix that lists the first (start) and last (end) row
for each factor in the sorted pattern matrix.}
}
\description{
faSort takes an unsorted factor pattern or structure matrix and returns a
sorted matrix with (possibly) reflected columns. Sorting is done such that
variables that load on a common factor are grouped together for ease of
interpretation.
}
\examples{

set.seed(123)
F <- matrix( c( .5,  0, 
                .6,  0,
                 0, .6,
                .6,  0,
                 0, .5,
                .7,  0,
                 0, .7,
                 0, .6), nrow = 8, ncol = 2, byrow=TRUE)

Rex1 <- F \%*\% t(F); diag(Rex1) <- 1

Items <- c("1. I am often tense.\n",
           "2. I feel anxious much of the time.\n",
           "3. I am a naturally curious individual.\n",
           "4. I have many fears.\n",
           "5. I read many books each year.\n",
           "6. My hands perspire easily.\n",
           "7. I have many interests.\n",
           "8. I enjoy learning new words.\n")

exampleOut <- fals(R = Rex1, nfactors = 2)

# Varimax rotation
Fload <- varimax(exampleOut$loadings)$loadings[]

# Add some row labels
rownames(Fload) <- paste0("V", 1:nrow(Fload))

cat("\nUnsorted fator loadings\n")
print(round( Fload, 2) )

# Sort items and reflect factors
out1 <- faSort(fmat = Fload, 
               salient = .25, 
               reflect = TRUE)
               
FloadSorted <- out1$loadings

cat("\nSorted fator loadings\n")
print(round( FloadSorted, 2) )

# Print sorted items
cat("\n Items sorted by Factor\n")
cat("\n",Items[out1$sortOrder])
}
\seealso{
\code{\link{fals}}

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
\keyword{Statstics}
