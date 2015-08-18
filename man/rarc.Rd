\name{rarc}
\alias{rarc}
\title{Rotate Points on the Surface on an N-Dimensional Ellipsoid}
\description{
Rotate between two points on the surface on an n-dimensional ellipsoid. The hyper-ellipsoid is composed of all points, B, such that B' Rxx B = Rsq. Vector B contains standardized regression coefficients.
}
\usage{
rarc(Rxx, Rsq, b1, b2, Npoints)
}
\arguments{
 \item{Rxx}{Predictor correlation matrix.}
  \item{Rsq}{Model coefficient of determination.}
  \item{b1}{First point on ellipsoid. If b1 and b2 are scalars then choose scaled 
           eigenvectors v[b1] and v[b2] as the start and end vectors.}
  \item{b2}{Second point on ellipsoid. If b1 and b2 are scalars then choose scaled 
           eigenvectors v[b1] and v[b2] as the start and end vectors.}
  \item{Npoints}{Generate ``Npoints'' +1 OLS coefficient vectors between b1 and b2.}
}
\value{
 \item{b}{N+1 sets of OLS coefficient vectors between b1 and b2.} 
}
\references{
Waller, N. G. & Jones, J. A. (2011). Investigating the performance of alternate regression weights by studying all possible criteria in regression models with a fixed set of predictors. \emph{Psychometrika, 76}, 410-439.
}
\author{Niels Waller and Jeff Jones.}
\examples{

## Example
## GRE/GPA Data
##-------------------##
R <- Rxx <- matrix(c(1.00, .56, .77,
                      .56, 1.00, .73,
                      .77, .73, 1.00), 3, 3)
                  
## GPA validity correlations                 
rxy <- c(.39, .34, .38)
b <- solve(Rxx) \%*\% rxy
 
Rsq <- t(b) \%*\% Rxx \%*\% b
N <- 200       
                   
b <- rarc(Rxx = R, Rsq, b1 = 1, b2 = 3, Npoints = N) 
 
## compute validity vectors
r <- Rxx \%*\% b
N <- N + 1
Rsq.r <- Rsq.unit <- rep(0, N)

for(i in 1:N){
## eval performance of unit weights
  Rsq.unit[i] <- (t(sign(r[,i])) \%*\% r[,i])^2 /
 		           (t(sign(r[,i])) \%*\% R \%*\% sign(r[,i]))
 		           
## eval performance of correlation weights               
  Rsq.r[i] <- (t(r[,i]) \%*\% r[,i])^2 /(t(r[,i]) \%*\% R \%*\% r[,i])	
}
 
cat("\nAverage relative performance of unit weights across elliptical arc:",
 	    round(mean(Rsq.unit)/Rsq,3) )     
cat("\n\nAverage relative performance of r weights across elliptical arc:",
 	    round(mean(Rsq.r)/Rsq,3) ) 
 

plot(seq(0, 90, length = N), Rsq.r, typ = "l", 
          ylim = c(0, .20),
          xlim = c(0, 95),
          lwd = 3,
          ylab = expression(R^2),
          xlab = expression(paste("Degrees from ",b[1]," in the direction of ",b[2])),
          cex.lab = 1.25, lab = c(10, 5, 5))
 points(seq(0, 90, length = N), Rsq.unit, 
          type = "l", 
          lty = 2, lwd = 3)
 legend(x = 0,y = .12,
        legend = c("r weights", "unit weights"), 
        lty = c(1, 2),
        lwd = c(4, 3),
        cex = 1.5)               
}
\keyword{datagen}% at least one, from doc/KEYWORDS
