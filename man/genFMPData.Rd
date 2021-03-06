% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genFMPData.R
\name{genFMPData}
\alias{genFMPData}
\title{Generate item response data for a filtered monotonic polynomial IRT model}
\usage{
genFMPData(NSubj, bParams, theta = NULL, thetaMN = 0, thetaVar = 1, seed)
}
\arguments{
\item{NSubj}{the desired number of subject response vectors.}

\item{bParams}{a p(items)-by-9 matrix of polynomial coefficients and model
designations. Columns 1 - 8 hold the polynomial coefficients; column 9 holds
the value of \code{k}.}

\item{theta}{A user-supplied vector of latent trait scores. Default theta =
NULL.}

\item{thetaMN}{If theta = NULL genFMPdata will simulate random normal
deviates from a population with mean thetaMN and variance thetaVar.}

\item{thetaVar}{If theta = NULL genFMPData will simulate random normal
deviates from a population with mean thetaMN and variance thetaVar.}

\item{seed}{initial seed for the random number generator.}
}
\value{
\item{theta}{theta values used for data generation}
\item{data}{N(subject)-by-p(items) matrix of item response data.}
\item{seed}{Value of the random number seed.}
}
\description{
Generate item response data for the filtered polynomial IRT model.
}
\examples{


# The following code illustrates data generation for 
# an FMP of order 3 (i.e., 2k+1)

# data will be generated for 2000 examinees
NSubjects <- 2000


## Example item paramters, k=1 FMP 
b <- matrix(c(
    #b0    b1     b2    b3      b4   b5 b6 b7  k
  1.675, 1.974, -0.068, 0.053,  0,  0,  0,  0, 1,
  1.550, 1.805, -0.230, 0.032,  0,  0,  0,  0, 1,
  1.282, 1.063, -0.103, 0.003,  0,  0,  0,  0, 1,
  0.704, 1.376, -0.107, 0.040,  0,  0,  0,  0, 1,
  1.417, 1.413,  0.021, 0.000,  0,  0,  0,  0, 1,
 -0.008, 1.349, -0.195, 0.144,  0,  0,  0,  0, 1,
  0.512, 1.538, -0.089, 0.082,  0,  0,  0,  0, 1,
  0.122, 0.601, -0.082, 0.119,  0,  0,  0,  0, 1,
  1.801, 1.211,  0.015, 0.000,  0,  0,  0,  0, 1,
 -0.207, 1.191,  0.066, 0.033,  0,  0,  0,  0, 1,
 -0.215, 1.291, -0.087, 0.029,  0,  0,  0,  0, 1,
  0.259, 0.875,  0.177, 0.072,  0,  0,  0,  0, 1,
 -0.423, 0.942,  0.064, 0.094,  0,  0,  0,  0, 1,
  0.113, 0.795,  0.124, 0.110,  0,  0,  0,  0, 1,
  1.030, 1.525,  0.200, 0.076,  0,  0,  0,  0, 1,
  0.140, 1.209,  0.082, 0.148,  0,  0,  0,  0, 1,
  0.429, 1.480, -0.008, 0.061,  0,  0,  0,  0, 1,
  0.089, 0.785, -0.065, 0.018,  0,  0,  0,  0, 1,
 -0.516, 1.013,  0.016, 0.023,  0,  0,  0,  0, 1,
  0.143, 1.315, -0.011, 0.136,  0,  0,  0,  0, 1,
  0.347, 0.733, -0.121, 0.041,  0,  0,  0,  0, 1,
 -0.074, 0.869,  0.013, 0.026,  0,  0,  0,  0, 1,
  0.630, 1.484, -0.001, 0.000,  0,  0,  0,  0, 1), 
  nrow=23, ncol=9, byrow=TRUE)  

# generate data using the above item paramters
data<-genFMPData(NSubj = NSubjects, bParams=b, seed=345)$data

}
\author{
Niels Waller
}
\keyword{statistics}
