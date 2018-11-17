# scale vector to unit length
	vnorm <- function(x) as.matrix(x/c(sqrt(t(x) %*% x)))