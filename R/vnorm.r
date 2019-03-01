# scale vector to unit length


#' Norm a Vector to Unit Length
#' 
#' Norm a vector to unit length.
#' 
#' 
#' @param x An n by 1 vector.
#' @return the scaled (i.e., unit length) input vector
#' @author Niels Waller
#' @keywords statistics
#' @export
#' @examples
#' 
#'  x <- rnorm(5)
#'  v <- vnorm(x)
#'  print(v)
#' 
	vnorm <- function(x) as.matrix(x/c(sqrt(t(x) %*% x)))
