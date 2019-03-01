# compute the cosine between two vectors


#' Compute the Cosine Between Two Vectors
#' 
#' Compute the cosine between two vectors.
#' 
#' 
#' @param x A p x 1 vector.
#' @param y A p x 1 vector.
#' @return \item{Cosine between x and y}{}
#' @keywords Statistics
#' @export
#' @examples
#' 
#' x <- rnorm(5)
#' y <- rnorm(5)
#' vcos(x, y)
#' 
vcos <- function(x,y){
    if(length(x)!=length(y)) stop(" x and y have different lengths")
     vnorm <- function(x) as.matrix(x/as.numeric(sqrt(t(x) %*% x)))
     t(vnorm(x)) %*% vnorm(y)
}
