#' Calculate Univariate Skewness for a Vector or Matrix
#' 
#' Calculate univariate skewness for vector or matrix (algorithm G1 in Joanes &
#' Gill, 1998).
#' 
#' 
#' @param x Either a vector or matrix of numeric values.
#' @return \item{Skewness for each column in x.}{}
#' @author Niels Waller
#' @seealso \code{\link{kurt}}
#' @references Joanes, D. N. & Gill, C. A. (1998). Comparing measures of sample
#' skewness and kurtosis. \emph{The Statistician, 47}, 183-189.
#' @keywords Statistics
#' @export
#' @examples
#' 
#' x <- matrix(rnorm(1000), 100, 10)
#' skew(x)
#' 
skew<-function(x){
                     
          sk <- function(xx){
             n <- length(xx)
             mn <- mean(xx)
             dif.x <- xx - mn
             m2 <- sum(dif.x^2)/n
             m3 <- sum(dif.x^3)/n
             m4 <- sum(dif.x^4)/n
             b1 <- m3/(m2^(3/2))
             g1 <- (b1 * sqrt(n * (n - 1)))/(n - 2)
             g1
           }
          
          if(ncol(x)==1 || is.null(dim(x)))
             return(sk(x))
          else
           return(apply(x,2,sk))              
            
}
