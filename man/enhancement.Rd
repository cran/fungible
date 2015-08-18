\name{enhancement}
\alias{enhancement}
\title{Find OLS Regression Coefficients that Exhibit Enhancement}
\description{
Find OLS regression coefficients that exhibit a specified degree of enhancement. 
}
\usage{ 
enhancement(R, br, rr)
}
\arguments{
 \item{R}{Predictor correlation matrix.}
 \item{br}{Model R-squared = b' r. That is, br is the model coefficient of determination: b'Rb= Rsq = br}
 \item{rr}{Sum of squared predictor-criterion correlations (rxy). That is, rr = r'r = Sum(rxy^2)} 
}
\value{
 \item{b}{Vector of standardized regression coefficients.}
 \item{r}{Vector of predictor-criterion correlations.} 
}

\references{
Waller, N. G. (2011). The geometry of enhancement in multiple regression. \emph{Psychometrika, 76}, 634--649.
}
\author{Niels Waller}
\examples{
## Example: For a given predictor correlation  matrix (R) generate 
## regression coefficient vectors that produce enhancement (br - rr > 0)

## Predictor correlation matrix
R <- matrix(c( 1,  .5, .25,
              .5, 1,   .30,
              .25, .30, 1), 3, 3) 
 
## Model coefficient of determination
Rsq <- .60
 
output<-enhancement(R, br = Rsq, rr =.40) 
 
r <- output$r
b <- output$b
  
##Standardized regression coefficients
print(t(b)) 

##Predictor-criterion correlations
print(t(r)) 
 
##Coefficient of determinations (b'r)
print(t(b) \%*\% r)

##Sum of squared correlations (r'r)
print(t(r) \%*\% r)
}

\keyword{datagen}% at least one, from doc/KEYWORDS
