\name{HS9Var}
\alias{HS9Var}
\docType{data}
\title{9 Variables from the Holzinger and Swineford (1939) Dataset
}
\description{
Mental abilities data on seventh- and eighth-grade children  from the classic Holzinger and Swineford (1939) dataset.  
}
\usage{data("HS9Var")}
\format{
  A data frame with 301 observations on the following 15 variables.
  \describe{
    \item{id}{subject identifier}
    \item{sex}{gender}
    \item{ageyr}{age, year part}
    \item{agemo}{age, month part}
    \item{school}{school name (Pasteur or Grant-White)}
    \item{grade}{grade}
    \item{x1}{Visual perception}
    \item{x2}{Cubes}
    \item{x3}{Lozenges}
    \item{x4}{Paragraph comprehension}
    \item{x5}{Sentence completion}
    \item{x6}{Word meaning}
    \item{x7}{Speeded addition}
    \item{x8}{Speeded counting of dots}
    \item{x9}{Speeded discrimination straight and curved capitals}
  }
}
\source{
These data were retrieved from the lavaan package. The complete data for all 26 tests are available in the MBESS package. 
}
\references{
Holzinger, K., and Swineford, F. (1939). A study in factor analysis: The stability of a bifactor solution. Supplementary Educational Monograph, no. 48. Chicago: University of Chicago Press.

Joreskog, K. G. (1969). A general approach to confirmatory maximum likelihood factor analysis. \emph{Psychometrika, 34}, 183-202.
}
\examples{
data(HS9Var)
head(HS9Var)
}
\keyword{datasets}
