 ##############################################
 #  Version October 25, 2012
 # updated December 20, 2017
 #
 #############################################
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
 
 
