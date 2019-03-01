##--------------------------------------#
##  October 31, 2012
##  updated August 12, 2015
##  updated  Sept 25 2015
##  updated July 26, 2017
##  Niels Waller
##
##  Smooth a non-positive definite R matrix by the method of 
##  Knol and Ten Berge 1991 Eq (27)
##--------------------------------------#



#' Smooth a Non PD Correlation Matrix
#' 
#' A function for smoothing a non-positive definite correlation matrix by the
#' method of Knol and Berger (1991).
#' 
#' 
#' @param R A non-positive definite correlation matrix.
#' @param eps Small positive number to control the size of the non-scaled
#' smallest eigenvalue of the smoothed R matrix. Default = 1E8 *
#' .Machine$double.eps
#' @return \item{Rsmoothed}{A Smoothed (positive definite) correlation matrix.}
#' @author Niels Waller
#' @references Knol, D. L., and Berger, M. P. F., (1991). Empirical comparison
#' between factor analysis and multidimensional item response
#' models.\emph{Multivariate Behavioral Research, 26}, 457-477.
#' @keywords Statistics
#' @export
#' @examples
#' 
#' ## choose eigenvalues such that R is NPD
#' l <- c(3.0749126,  0.9328397,  0.5523868,  0.4408609, -0.0010000)
#' 
#' ## Generate NPD R
#' R <- genCorr(eigenval = l, seed = 123)
#' print(eigen(R)$values)
#' 
#' #> [1]  3.0749126  0.9328397  0.5523868  0.4408609 -0.0010000
#' 
#' ## Smooth R
#' Rsm<-corSmooth(R, eps = 1E8 * .Machine$double.eps)
#' print(eigen(Rsm)$values)
#' 
#' #> [1] 3.074184e+00 9.326669e-01 5.523345e-01 4.408146e-01 2.219607e-08
#' 
corSmooth <- function(R, eps = 1E8 * .Machine$double.eps){
  
     KDK <- eigen(R)
     Dplus <- D <- KDK$values
     Dplus[D <= 0]<- eps
     Dplus <- diag(Dplus)
     K <- KDK$vectors
    
     Rtmp <- K %*% Dplus %*% t(K)
     invSqrtDiagRtmp <- diag(1/sqrt(diag(Rtmp)))
     R<- invSqrtDiagRtmp %*% Rtmp %*% invSqrtDiagRtmp
     # for perfect symmetry
     .5*( R + t(R) )
}

