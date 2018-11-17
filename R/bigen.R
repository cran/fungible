##########################################################
# This function generates binary data with user-defined  #
# thresholds                                             #
#                                                        #
# Arguments:                                             #
# data - Either a matrix of binary (0/1) indicators or a #
#       correlation matrix.                              #
#                                                        #
# n    - The desired sample size of the simulated data.  #
#                                                        #
# thresholds  - If data is a correlation matrix then     #
#               thresholds must be a vector of threshold #
#               cut points.                              #
# Smooth      - (logical) smooth = TRUE will smooth the  #
#               tetrachoric correltion matrix            #
#                                                        #
# Output                                                 #
# data	-       Simulated binary data                    #
#  r	-         Input or calculated (tetrachoric)        #
#               correlation matrix                       #
#                                                        #
##########################################################



bigen <- function(data, n, thresholds = NULL, 
                  Smooth = FALSE, seed = NULL){

 if(!is.null(seed)) set.seed(seed)
  
 nr <- nrow(data)
 nitems <- ncol(data)


##-----------------
## Function rmvnorm by F. Leisch
## a random number generator for the multivariate normal 
## distribution with mean equal to mean and covariance matrix sigma.
 rmvnorm<- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean))) 
 {
    if (nrow(sigma) != ncol(sigma)) {
        stop("sigma must be a square matrix")
    }
    if (length(mean) != nrow(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    sigsvd <- svd(sigma)
    retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
    retval <- sweep(retval, 2, mean, "+")
    retval
 } #END rmvrnorm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
## If binary data supplied, compute population tetrachoric 
## correlation matrix 
 if(nr > nitems){
  ## compute thresholds from supplied binary data
     thresholds <- qnorm(apply(data, 2, mean))
     bidat <- matrix(0, nrow = n, ncol = nitems)
     r <- tetcor(X = data, stderror = F, Smooth = Smooth)$r
  ## Generate MVN data     
     ghost <- rmvnorm(n, mean = rep(0, nitems), sigma = r)
  ## dichotomize data at thresholds   
     for(i in 1:nitems){
         bidat[ghost[,i] <= thresholds[i], i] <- 1
     }
    
    result <- list(data = bidat, r=r) 
 }# END if(nr > nitems)

## If population correlation matrix supplied
 if(nr == nitems){

     if(is.null(thresholds))stop("thresholds must be supplied with r matrix input")

     bidat <- matrix(0,nrow=n,ncol=nitems) 
     
 ##  if R = I
     if(sum(data) == nitems){
       ghost <- matrix(rnorm(n*nitems), n, nitems)
     }  
     else {    # R != I
       ghost <- rmvnorm(n, mean = rep(0,nitems), sigma = data)
     }
  ## dichotomize data at thresholds   
     for(i in 1:nitems){
         bidat[ghost[, i] <= thresholds[i],i] <- 1
     }   
    result<-list(data = bidat, r = data) 
  }
 result
 } #END fnc bigen
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##


# N <- 500
# R <- matrix(.5,5,5)
# diag(R) <- 1
# 
# # desired difficulties (i.e., means) of 
# # the dichotomized scores
# difficulties <- c(.2, .3, .4, .5, .6)
# 
# # cut the observed scores at these thresholds
# # to approximate the above difficulties
# thresholds <- qnorm(difficulties)
# 
# 
# ## Now use 'bigen' to generate binary data matrix with 
# ## same correlations as in Binary
# 
# out <- bigen(data = R, n = N, thresholds = thresholds)
# 
# tetcor(out$data, stderror = F, PRINT = F)$r
# apply(out$data, 2, mean)

