#' Smooth a Non PD Correlation Matrix using the Knol-Berger algorithm
#' 
#' A function for smoothing a non-positive definite correlation matrix by the
#' method of Knol and Berger (1991).
#' 
#' 
#' @param R A non-positive definite correlation matrix.
#' @param eps Small positive number to control the size of the non-scaled
#' smallest eigenvalue of the smoothed R matrix. Default = 1E8 *
#' .Machine$double.eps
#' @return \item{RKB}{A Smoothed (positive definite) correlation matrix.}
#' \item{eps}{Small positive number to control the size of the non-scaled
#' smallest eigenvalue of the smoothed R matrix.}
#' @author Niels Waller
#' @references Knol, D. L., & Berger, M. P. F., (1991). Empirical comparison
#' between factor analysis and multidimensional item response
#' models.\emph{Multivariate Behavioral Research, 26}, 457-477.
#' @keywords Statistics
#' @export
#' @examples
#' 
#' data(BadRLG)
#' 
#' ## RKB = smoothed R
#' RKB<-smoothKB(R=BadRLG, eps = 1E8 * .Machine$double.eps)$RKB
#' print(eigen(RKB)$values)
#' 
#' 
smoothKB <- function(R, eps = 1E8 * .Machine$double.eps){
  ##--------------------------------------#
  ##  October 31, 2012
  ##  updated August 12, 2015
  ##  updated  Sept 25 2015
  ##  updated May 18, 2019
  ##  Niels Waller
  ##
  ##  Smooth a non-positive definite R matrix by the method of 
  ##  Knol and Berger 1991 Eq (27)
  ##--------------------------------------#
     Nvar <- nrow(R)
     KDK <- eigen(R)
     Dplus <- D <- KDK$values
     Dplus[D <= 0]<- eps
     Dplus <- diag(Dplus)
     K <- KDK$vectors
    
     Rtmp <- K %*% Dplus %*% t(K)
     invSqrtDiagRtmp <- diag(1/sqrt(diag(Rtmp)))
     RKB <- invSqrtDiagRtmp %*% Rtmp %*% invSqrtDiagRtmp
     
     # insure perfect symmetry
     RKB <- .5 * (RKB + t(RKB))
      
     colnames(RKB) <- rownames(RKB) <- colnames(R)

	 list(RKB = RKB, eps = eps)
}

