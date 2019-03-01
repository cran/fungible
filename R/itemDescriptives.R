 ##############################################
 #  Version October 25, 2012
 # updated December 20, 2017
 #
 #############################################


#' Compute basic descriptives for binary-item analysis
#' 
#' Compute basic descriptives for binary item analysis
#' 
#' 
#' @param X a matrix of binary (0/1) item responses.
#' @param digits number of digits to print.
#' @return \item{alpha}{Coefficient alpha for the total scale.}
#' \item{means}{item means.} \item{standard deviations}{item standard
#' deviations.} \item{pt. biserial correlations}{corrected item-total point
#' biserial correlations.} \item{biserial correlations}{corrected item-total
#' point biserial correlations.} \item{corrected.alpha}{corrected (leave item
#' out) alpha coefficients.}
#' @author Niels Waller
#' @keywords statistics
#' @export
#' @examples
#' 
#' 	## Example 1: generating binary data to match
#' 	## an existing binary data matrix
#' 	##
#' 	## Generate correlated scores using factor 
#' 	## analysis model
#' 	## X <- Z *L' + U*D 
#' 	## Z is a vector of factor scores
#' 	## L is a factor loading matrix
#' 	## U is a matrix of unique factor scores
#' 	## D is a scaling matrix for U
#' 
#' 	Nsubj <- 2000
#' 	L <- matrix( rep(.707,5), nrow = 5, ncol = 1)
#' 	Z <-as.matrix(rnorm(Nsubj))
#' 	U <-matrix(rnorm(Nsubj * 5),nrow = Nsubj, ncol = 5)
#' 	tmp <-  sqrt(1 - L^2) 
#' 	D<-matrix(0, 5, 5)
#' 	diag(D) <- tmp
#' 	X <- Z %*% t(L) + U%*%D
#' 
#' 	cat("\nCorrelation of continuous scores\n")
#' 	print(round(cor(X),3))
#' 
#' 	thresholds <- c(.2,.3,.4,.5,.6)
#' 
#' 	Binary<-matrix(0,Nsubj,5)
#' 	for(i in 1:5){
#' 	  Binary[X[,i]<=thresholds[i],i]<-1
#' 	}   
#' 
#' 	cat("\nCorrelation of Binary scores\n")
#' 	print(round(cor(Binary),3))
#' 
#' 	## Now use 'bigen' to generate binary data matrix with 
#' 	## same correlations as in Binary
#' 
#' 	z <- bigen(data = Binary, n = 5000)
#' 
#' 	cat("\n\nnames in returned object\n")
#' 	print(names(z))
#' 
#' 	cat("\nCorrelation of Simulated binary scores\n")
#' 	print(round( cor(z$data), 3))
#' 
#' 
#' 	cat("Observed thresholds of simulated data:\n")
#' 	cat( apply(z$data, 2, mean) )
#' 	
#' 	itemDescriptives(z$data)
#' 
 itemDescriptives <- function(X, digits = 3){
       total  <- apply(X, 1, sum)
       Nitems <- ncol(X)
       Nsubj  <- nrow(X)

       mns <- apply(X, 2, mean)
       sds <- apply(X, 2, sd)
	   
	     Alpha <- function(covar){
	       if(dim(covar)[1] !=dim(covar)[2]) stop("Input matrix must be Square Covar")
	       k <- nrow(covar)
	       (k/(k - 1)) * (1 - sum(diag(covar))/sum(covar))
	     }

       ##Compute corrected item total pt.biserial and correlcted alpha
       corrected.pt.biserial <- corrected.alpha <- rep( 0, ncol(X) )
       for(i in 1:Nitems){
         corrected.pt.biserial[i] <- cor( X[,i], (total-X[,i]) )
         corrected.alpha[i] <- Alpha( var(X[,-i]) ) 
       }


       ##Compute corrected  item total biseral
       corrected.biserial <- rep(0, ncol(X))
       for(i in 1:Nitems){
          total.corrected <- total - X[,i]
          mu.plus <- mean( total.corrected[X[,i] == 1] )
          mu.x <- mean(total.corrected)
          sd.x <- sd(total.corrected)
          p <- mns[i]
          y <- dnorm(qnorm(p))
          corrected.biserial[i] <- (p/y) * (mu.plus - mu.x)/sd.x
          }

        cat("\n\nReliability: alpha = ",round(Alpha(var(X)),2),"\n\n")
        results <- data.frame(Means = round(mns,digits),
                   SDs = round(sds,digits),
                   corrected.pt.biserial=round(corrected.pt.biserial,digits),
                   corrected.biserial=round(corrected.biserial,digits),
                   corrected.alpha = round(corrected.alpha, digits))
        print(results)
        invisible(results)

 }
 
 
