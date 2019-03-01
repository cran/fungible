##---------------------------------------------#
#  Generate sample correlation matrices
#  see Browne 1968
#  Browne, M. (1968). A comparison of factor analytic techniques. 
#  Psychometrika, 33(3):267-334.
#
# Output
#    cor.sample    a sampled correlation matrix
#    cov.sample    a sampled covariance matrix
#----------------------------------------------#



#' Sample Correlation Matrices from a Population Correlation Matrix
#' 
#' Sample correlation (covariance) matrices from a population correlation
#' matrix (see Browne, 1968; Kshirsagar, 1959)
#' 
#' 
#' @param R A population correlation matrix.
#' @param n Sample correlation (covariance) matrices will be generated assuming
#' a sample size of n.
#' @return \item{cor.sample}{Sample correlation matrix.}
#' \item{cov.sample}{Sample covariance matrix.}
#' @author Niels Waller
#' @references Browne, M. (1968). A comparison of factor analytic techniques.
#' \emph{Psychometrika, 33(3)}, 267-334.
#' 
#' Kshirsagar, A. (1959). Bartlett decomposition and Wishart distribution.
#' \emph{The Annals of Mathematical Statistics, 30(1)}, 239-241.
#' @keywords datagen
#' @examples
#' 
#' R <- matrix(c(1, .5, .5, 1), 2, 2)
#' # generate a sample correlation from pop R with n = 25
#' out <- corSample(R, n = 25)
#' out$cor.sample
#' out$cov.sample
#' @importFrom stats rnorm rchisq 
#' @export

corSample<- function(R,n){
  Nvar<-ncol(R)
	Tmat<-matrix(0,Nvar,Nvar)
	Tmat[lower.tri(Tmat)]<-rnorm(n=(Nvar*(Nvar-1)/2))
	for(i in 1:Nvar){
# Note that Browne 1968 contains a typo for the df -- see the following for (n-i+1)
# Kshirsagar, A. (1959). Bartlett decomposition and wishart distribution. The Annals of Mathematical Statistics, 
#    30(1)239-241.
	 Tmat[i,i]<-sqrt(rchisq(n=1,df=(n-i+1)))
	}

	H<- Tmat %*% t(Tmat)
	
	Omega <-t(chol(R))
	A <-  Omega %*% H %*% t(Omega)
	S <- (1/n) * A
	Dmat<-diag(1/sqrt(diag(A)))
	R.samp<-Dmat%*%A%*%Dmat
	R.samp <- .5*(R.samp + t(R.samp))
	list(cor.sample=R.samp, cov.sample=S)
}


