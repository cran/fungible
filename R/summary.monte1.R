#' Summary Method for an Object of Class Monte1
#' 
#' summary method for class "monte1"
#' 
#' 
#' @aliases summary.monte1 print.summary.monte1
#' @param object An object of class \code{monte1}, usually, a result of a call
#' to \code{monte1}.
#' @param digits Number of significant digits to print in final results.
#' @param \dots Additional argument affecting the summary produced.
#' @return Various descriptive statistics will be computed including"
#' \enumerate{
#'   \item{Expected correlation matrix.}
#'   \item{Observed correlation matrix.}
#'   \item{Expected indicator skewness values.}
#'   \item{Observed indicator skewness values.}
#'   \item{Expected indicator kurtosis values.}
#'   \item{Observed indicator kurtosis values.}
#' }
#' @export
#' @examples
#' 
#' ## Generate dimensional data for 4 variables. 
#' ## All correlations = .60; all variable
#' ## skewness = 1.75; 
#' ## all variable kurtosis = 3.75
#' 
#' cormat <- matrix(.60, 4, 4)
#' diag(cormat) <- 1
#' 
#' nontaxon.dat <- monte1(seed = 123, nsub = 100000, nvar = 4, skewvec = rep(1.75, 4),
#'                  kurtvec = rep(3.75, 4), cormat = cormat)
#' 
#' summary(nontaxon.dat)
#' 
summary.monte1<-function(object, digits=3, ...){
z<-object
nvar<-z$nvar


obs.skew.matrix<-obs.kurt.matrix<-matrix(99,1,nvar)                                  
                

  obs.cormat  <- cor(z$data)
  obs.skewvec <- apply(z$data, 2, skew)
  obs.kurtvec <- apply(z$data, 2, kurt)


ans <- z["call"]


cat("\nCall monte1:","\n")
print(z$call)

cat("\nSeed = ",z$seed,"\n")



cat("\n\n\nExpected correlation matrix:","\n")
print(round(z$cormat,digits))

cat("\n\nObserved correlation matrix:","\n")
print(round(obs.cormat,digits))


cat("\nExpected indicator skewness:")
    print(round(z$skewvec, digits))
   
cat("\nObserved indicator skewness:","\n")
print(round(obs.skewvec,digits))

cat("\n\nExpected indicator kurtosis:")
    print(round(z$kurtvec,digits))

cat("\nObserved indicator kurtosis:","\n")
 print(round(obs.kurtvec,digits))


 
 ans$seed<-z$seed
 ans$nvar <- nvar
 ans$cormat  <- z$cormat
 ans$obs.cormat<-obs.cormat
 ans$skewvec <- z$skewvec
 ans$obs.skewvec <-obs.skewvec
 ans$kurtvec <- z$kurtvec
 ans$obs.kurtvec <- obs.kurtvec
 
 
 class(ans)<-"summary.monte1"
invisible(ans)
}
