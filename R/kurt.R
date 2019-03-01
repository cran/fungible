#######################################################
# algorithm G2 in D. N. Joanes and C. A. Gill (1998), 
# Comparing measures of sample skewness and kurtosis. The Statistician, 
# 47, 183-189.



#' Calculate Univariate Kurtosis for a Vector or Matrix
#' 
#' Calculate univariate kurtosis for a vector or matrix (algorithm G2 in Joanes
#' & Gill, 1998).
#' 
#' 
#' @param x Either a vector or matrix of numeric values.
#' @return \item{Kurtosis for each column in x.}{}
#' @author Niels Waller
#' @seealso \code{\link{skew}}
#' @references Joanes, D. N. & Gill, C. A. (1998). Comparing measures of sample
#' skewness and kurtosis. \emph{The Statistician, 47}, 183-189.
#' @keywords Statistics
#' @export
#' @examples
#' 
#' x <- matrix(rnorm(1000), 100, 10)
#' print(kurt(x))
#' 
kurt<-function(x){
 
       kt<-function(xx){
            n <- length(xx)
            mn <- mean(xx)
            dif.x <- xx - mn
            m2 <- sum(dif.x^2)/n
            m4 <- sum(dif.x^4)/n
            b2 <- (m4/m2^2)       
            g2 <- ((n + 1) * (n - 1) * (b2 - (3 * (n - 1))/(n + 1)))/((n -2) * (n - 3))
            g2
       }     
                
           if(ncol(x)==1 || is.null(dim(x)))
             return(kt(x))
          else
           return(apply(x,2,kt))                                                      
 }         
