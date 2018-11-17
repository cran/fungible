\name{bigen}
\alias{bigen}
\title{Generate Correlated Binary Data}
\description{
  Function for generating binary data with population thresholds.
}
\usage{
bigen(data, n, thresholds = NULL, Smooth = FALSE, seed = NULL)
}
\arguments{
 \item{data}{
  Either a matrix of binary (0/1) indicators or a correlation matrix. 
 }
 \item{n}{
  The desired sample size of the simulated data. 
 }
 \item{thresholds}{
 If data is a correlation matrix, thresholds must be a vector of threshold cut points.}
 \item{Smooth}{
 (logical) Smooth = TRUE will smooth the tetrachoric correltion matrix.}
 \item{seed}{Default = FALSE. Optional seed for random number generator.}
}
\value{
  \item{data}{Simulated binary data}
  \item{r}{Input or calculated (tetrachoric) correlation matrix}
}
\author{Niels G Waller}
\examples{
## Example: generating binary data to match
## an existing binary data matrix
##
## Generate correlated scores using factor 
## analysis model
## X <- Z *L' + U*D 
## Z is a vector of factor scores
## L is a factor loading matrix
## U is a matrix of unique factor scores
## D is a scaling matrix for U

N <- 5000

# Generate data from a single factor model
# factor patter matrix
L <- matrix( rep(.707, 5), nrow = 5, ncol = 1)

# common factor scores
Z <- as.matrix(rnorm(N))

# unique factor scores
U <- matrix(rnorm(N *5), nrow = N, ncol = 5)
D <- diag(as.vector(sqrt(1 - L^2)))

# observed scores
X <- Z \%*\% t(L) + U \%*\% D

cat("\nCorrelation of continuous scores\n")
print(round(cor(X),3))

# desired difficulties (i.e., means) of 
# the dichotomized scores
difficulties <- c(.2, .3, .4, .5, .6)

# cut the observed scores at these thresholds
# to approximate the above difficulties
thresholds <- qnorm(difficulties)

Binary <- matrix(0, N, ncol(X))
for(i in 1:ncol(X)){
  Binary[X[,i] <= thresholds[i],i] <- 1
}   

cat("\nCorrelation of Binary scores\n")
print(round(cor(Binary), 3))

## Now use 'bigen' to generate binary data matrix with 
## same correlations as in Binary

z <- bigen(data = Binary, n = N)

cat("\n\nnames in returned object\n")
print(names(z))

cat("\nCorrelation of Simulated binary scores\n")
print(round(cor(z$data), 3))


cat("Observed thresholds of simulated data:\n")
cat(apply(z$data, 2, mean))
}
\keyword{datagen}
