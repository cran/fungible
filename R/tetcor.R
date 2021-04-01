#' Compute ML Tetrachoric Correlations
#' 
#' Compute ML tetrachoric correlations with optional bias correction and
#' smoothing.
#' 
#' 
#' @param X Either a matrix or vector of (0/1) binary data.
#' @param y An optional(if X is a matrix) vector of (0/1) binary data.
#' @param BiasCorrect A logical that determines whether bias correction (Brown
#' & Benedetti, 1977) is performed. Default = TRUE.
#' @param stderror A logical that determines whether standard errors are
#' calulated.  Default = FALSE.
#' @param Smooth A logical which determines whether the tetrachoric correlation
#' matrix should be smoothed.  A smoothed matrix is always positive definite.
#' @param max.iter Maximum number of iterations. Default = 50.
#' @param PRINT A logical that determines whether to print progress updates
#' during calculations. Default = TRUE
#' @return If stderror = FALSE, \code{tetcor} returns a matrix of tetrachoric
#' correlations. If \code{stderror} = TRUE then \code{tetcor} returns a list
#' the first component of which is a matrix of tetrachoric correlations and the
#' second component is a matrix of standard errors (see Hamdan, 1970).
#' 
#' @return 
#'  \item{r}{The tetrachoric correlation matrix}.
#'  \item{se}{A matrix of standard errors.}
#'  \item{convergence}{(logical) The convergence status of the algorithm. A value of 
#'    TRUE denotes that the algorithm converged. A value of FALSE denotes that the 
#'    algorithm did not converge and the returned correlations 
#'    are Pearson product moments.}
#'  \item{Warnings}{A list of warnings.}
#' @author Niels Waller
#' @references Brown, M. B. & Benedetti, J. K. (1977). On the mean and variance
#' of the tetrachoric correlation coefficient. \emph{Psychometrika, 42},
#' 347--355.
#' 
#' Divgi, D. R. (1979) Calculation of the tetrachoric correlation coefficient.
#' \emph{Psychometrika, 44}, 169-172.
#' 
#' Hamdan, M. A. (1970). The equivalence of tetrachoric and maximum likelihood
#' estimates of rho in 2 by 2 tables. \emph{Biometrika, 57}, 212-215.
#' @keywords Statistics
#' @import mvtnorm
#' @export
#' @examples
#' 
#' ## generate bivariate normal data
#' library(MASS)
#' set.seed(123)
#' rho <- .85
#' xy <- mvrnorm(100000, mu = c(0,0), Sigma = matrix(c(1, rho, rho, 1), ncol = 2))
#' 
#' # dichotomize at difficulty values
#' p1 <- .7
#' p2 <- .1
#' xy[,1] <- xy[,1] < qnorm(p1) 
#' xy[,2] <- xy[,2] < qnorm(p2)
#' 
#' print( apply(xy,2,mean), digits = 2)
#' #[1] 0.700 0.099
#' 
#' tetcor(X = xy, BiasCorrect = TRUE, 
#'        stderror = TRUE, Smooth = TRUE, max.iter = 5000)
#' 
#' # $r
#' # [,1]      [,2]
#' # [1,] 1.0000000 0.8552535
#' # [2,] 0.8552535 1.0000000
#' # 
#' # $se
#' # [,1]           [,2]
#' # [1,] NA         0.01458171
#' # [2,] 0.01458171 NA
#' # 
#' # $Warnings
#' # list()
#' 
#' 
#' 
tetcor<-function (X, y = NULL, BiasCorrect = TRUE, stderror = FALSE, 
                  Smooth = TRUE, max.iter = 5000, PRINT = TRUE) 
{
# add: cparam=NULL to argument list
#------------------------------------------------------------#
# Function: tetcor                                           #
#                                                            #
# December 16, 2017                                          #
#    Improved fnc speed via FastTable                        #
#                                                            #
#  April 27, 2016                                            #
# fixed se's for cases in which two cells are empty          #
# May 1 2016 define cell counts as.numeric to avoid          #
#  integer overflow                                          #
#                                                            #
#  warning messages updated Sept 2 2004                      # 
#                                                            #
# Niels Waller, June 9, 2018                                 #
# Requires library:  mvtnorm from CRAN                       #
#                                                            #
# Compute ML estimate of tetrachoric r                       #
#                                                            #
# ARGUMENTS:                                                 #
#  X :  a matrix of binary responses coded 0/1               #
#  y : if x is a vector then y must be a vector of 0/1       #
#     responses                                              #
#  BiasCorrect : a logical that determines whether to apply  #
#  the Brown & Benedetti bias correction                     #
#  stderror : a logical that determines whether to compute   #
#             the s.e. of the estimate                       # 
#                                                            #
#  Smooth :  a logical that determines whether to smooth the #
#            correlation matrix                              #
#                                                            #
#  max.iter :                                                #
#                                                            #
#                                                            #
# VALUE                                                      #
#   r :   matrix of tetrachoric correlations                 #
#   se : matrix of standard error of r is requested          #
#   Warnings : a list of correlations that did not converge  #
#                                                            #
#                                                            #
#                                                            #
# Initial estimate of r calculated by Divgi's method         #
#                                                            #
# Divgi, D. R. (1979) Calculation of the tetrachoric         #
#  correlation coefficient. Psychometrika, 44, 169--172.     #
#                                                            #
# Bias correction:                                           #
# Brown, M. B. & Benedetti, J. K. (1977). On the mean and    #
#  variance of the tetrachoric correlation coefficient.      #
#  Psychometrika, 42, 347--355.                              #
#                                                            #
# S.E. of r_{tet}                                            #
# Hamdan, M. A. (1970). The equivalence of tetrachoric and   #
#  maximum likelihood estimates of $\rho$ in $2 \times 2$    #
#  tables. Biometrika, 57, 212--215.                         #
#----------------------------------------------------------- #

# initialize warning flag
 warning.msg <- list()
 
 x <- X

# Update: September 8, 2014
# use 'cuhre' for bivariate integral 
# library(R2Cuba)

  if (is.data.frame(x)) 
      x <- as.matrix(x) 
  if (is.data.frame(y)) 
       y <- as.matrix(y)
  if (!is.matrix(x) && is.null(y)) 
        stop("supply both x and y or a matrix-like x")
  if(!is.null(y)) 
     x<-cbind(x,y)
  if(sum(is.na(x)) > 0)
        stop("Missing values are not allowed")
     
  nitems <- ncol(x)
  N <- nrow(x)
  
  #--Do any items have zero variance?
  zeroVar <- apply(x, 2, var) == 0
  
  if (sum(zeroVar) > 0) {
    labs <- paste0("Item", 1:nitems)
    stop(paste0("\nSERIOUS ERROR \nThe following item(s) has variance = 0: ", 
                labs[zeroVar], "\n"))
  }
  
  #--Initialize matrices
  # an element of seFlagMatrix will be set to 99
  # if the aymptotic standard error can not 
  #  be computed due to low cell counts
  r <- matrix(0, nitems, nitems)
  se <- 99
  if(stderror == TRUE) se <- seFlagMatrix <- r
  

 warning.k <- 0  


## Function definitions
#-----compute bivariate normal density----------#
  bvn <- function(z){
     1/(2*pi*sqrt(1-rhat^2)) * exp(-(z[1]^2+z[2]^2-2*rhat*z[1]*z[2])
                                /(2*(1-rhat^2)) ) 
  }

  Fasttable <- function(x, y) {
    tt <- sum(x * y)
    t1 <- sum(x)
    t2 <- sum(y)
    matrix(c(length(x)-t1-t2+tt, t1-tt, t2-tt, tt), 2, 2)
  }
  
  

########################  
#--MAIN loop begins here
########################
 for(iter.row in 2:nitems){
    for(iter.col in 1:(iter.row-1)){
   
    if(PRINT){   
     cat(c("Working on correlation: ",iter.row, iter.col, "\n"))
    }   
          
    v1 <- x[,iter.row]
    v2 <- x[,iter.col]

    # we assume that data are 0/1 so this no longer needed
    #v1 <- factor(v1,levels=c("0","1"))
    #v2 <- factor(v2,levels=c("0","1"))
    
    abcd <- Fasttable(v1,v2)  
  
    # use as.numeric to avoid integer overflow when N is very large
    cella <- as.numeric(abcd[2,2])
    cellb <- as.numeric(abcd[2,1])
    cellc <- as.numeric(abcd[1,2])
    celld <- as.numeric(abcd[1,1])
    
#----Brown and Benedetti bias correction----------#

  if(BiasCorrect == TRUE){
    # see Brown and Benedetti page 353
    # when two cells are zero (either diagonal or off diagonal) assign
    # r of plus or minus 1
       if(cella == 0 & celld == 0){
         r[iter.row,iter.col] <- -1
         if(stderror == TRUE){
           seFlagMatrix[iter.row, iter.col] <-  seFlagMatrix[iter.col, iter.row] <- 99
         }
         next
       }
       if(cellc==0 & cellb==0){
         r[iter.row, iter.col] <- 1
         if(stderror == TRUE){
           seFlagMatrix[iter.row, iter.col] <- seFlagMatrix[iter.col, iter.row] <- 99
         }
          next
       }
       if(cella == 0 & cellb > 0 & cellc > 0 & celld > 0){   # only cell a is 0.00
         cella<-.5;  cellb<- cellb - .5;  cellc<- cellc - .5;  celld<- celld +.5      
       }
       if(cellb == 0 & cella > 0 & cellc > 0 & celld > 0){   # only cell b is 0.00
         cellb<-cellb+.5; cella<-cella - .5; cellc<-cellc+.5; celld<-celld-.5                
       }
      if(cellc == 0 & cellb > 0 & cella > 0 & celld > 0){   # only cell c is 0.00
         cellc<-cellc+.5; celld<-celld-.5; cella<-cella-.5; cellb<-cellb+.5
      } 
      if(celld == 0 & cellb > 0 & cellc > 0 & cella > 0){   # only cell d is 0.00 
         celld <- .5; cellb <- cellb - .5; cellc <- cellc - .5; cella <- cella + .5
      }
   }#-------------end of bias correction----------------------------#
 
    p1 <- (cella + cellb)/N
    p2 <- (cellc + cella)/N
        
    p11 <- cella/N
    p00 <- celld/N 
    p01 <- cellc/N
    p10 <- cellb/N   
    

  #-------Correction for non-zero lower asymptotes
  # Currently not working
  # Carroll
  # if(!is.null(cparam)){
  #   gi <- cparam[iter.row]
  #   wi <- 1-gi
  #   gj <- cparam[iter.col]
  #   wj <- 1-gj 
  # 
  #   t00 <- p00/(wi * wj)
  #   
  #   t01 <- (wj*p01 * gj*p00)/(wi*wj)
  #   if(t01 < 0) t01 <-0
  #   
  #   t10 <- (wi*p10 - gi*p00)/(wi*wj)
  #   if(t10 < 0) t10<-0
  #   
  #   t11 <- 1 - t00 - t01 - t10 
  #   
  #   p11 <- t11
  #   p1 <- 1 - ( (1-p1)/(1-gi) )     
  #   p2 <- 1 - ( (1-p2)/(1-gj) )  
  #   
  #   cella <- t11 * N
  #   cellb <- t10 * N
  #   cellc <- t01 * N
  #   celld <- t00 * N
  #     
  #}  # end correction for guessing  
        
    h.x <- qnorm(p1)
    k.y <- qnorm(p2)
       
    h.sign <- sign(h.x)
    k.sign <- sign(k.y)
    
    hstar <- max(abs(h.x), abs(k.y))
    kstar <- min(abs(h.x), abs(k.y))
    
    h2 <- hstar^2
    k2 <- kstar^2

   Aa <- .5/(1 + (h2+k2) * (.12454 - .27102  * ( 1 - hstar/sqrt(h2+k2) ) ) )
   Bb <- .5/(1 + (h2+k2) * (.82281 - 1.03514 * ( kstar/sqrt(h2+k2) ) ) )
   Cc <- .07557 * hstar + (hstar-kstar)^2 * (.51141/(hstar + 2.05793) -.07557/hstar)
   R  <- (cella * celld)/(cellb * cellc) 
   
   Dd <- h.sign * k.sign * kstar * (.79289 + 4.28981/(1+3.30231 * hstar))
   alpha <- Aa + Bb * (-1 + 1/(1+Cc * (log10(R) - Dd)^2))
    
          
   rhat<-cos(pi/(1 + R^alpha)) # rhat is start value
   if(is.na(rhat > 0)) rhat <- 0  # if x and y have means of exactly .5 rhat is NaN
     
   if( abs(rhat) == 1 ) {        # integration breaks down if |rhat| = 1
      rhat <-.95 * rhat
   }
 
#-----End Divgi's method for start values -----------------#

#------------------------
# Divgi sends the readers to Pearson 1900 for 
# the formula for this derivative. This form is
## reported in Kirk 1973, p. 261 Psychometrika
   h.xSq <- h.x^2
   k.ySq <- k.y^2
 
#---1st derivative of likelihood (a/N) wrt r------------#
 dLdr<-function(r){
      1/(2 * pi * sqrt(1-r^2)) * exp( - (h.xSq + k.ySq - 2 * h.x * k.y * r)/(2 * (1-r^2)) )
 }


 eps <- 99
 iterations <- 0
 
#-----Newton Raphson Loop to improve initial estimate
 
 converged <- TRUE

 while( eps > .00001){
#----integrate over bivariate normal surface to estimate (d/N)


### Old code uses 'cuhre' from the R2Cuba in place of adapt 
# adaptOut <- R2Cuba::cuhre(ndim = 2, ncomp = 1, integrand = bvn,
#          lower = c(-6, -6), upper = c(h.x, k.y), flags = list(verbose = 0),
#          rel.tol = .000001, abs.tol = 0, min.eval = 0, max.eval=max.iter)
# Now using mvtnorm
   
    adaptOut <- mvtnorm::pmvnorm(lower = c(-8, -8), 
                                 upper = c(h.x, k.y),
                                 sigma = matrix(c(1, rhat, rhat, 1),2,2))[1] 
                             
    
     ep11 <- adaptOut

     firstderiv <- dLdr(rhat)
     
     if( firstderiv <= 0.02 ) firstderiv <- firstderiv *  5.5 #pull deriv back if too low
     rnew <- rhat - (ep11-p11)/firstderiv
     
     if( abs(rnew) > 1 ) {
         rnew <- .99
         rhat <- 99 # do not stop during this iteration
     }
   
  
    eps <- abs(rnew - rhat)
    iterations <- iterations + 1 
  
    rhat <- rnew
    if(iterations > max.iter){
      converged <- FALSE
       warning.k <- warning.k + 1
       warn.msg <- paste("WARNING: Correlation ", iter.row, ", ", iter.col, 
                      " failed to converge!", sep="")
      # cat(warn.msg,"\n")
       warning.msg[warning.k] <- warn.msg
       rnew <- 51; 
       ## replace with Pearson
       rhat <- cor(x[,c(iter.row,iter.col)])[1,2]
       eps <- .00000001
      }   
    
 }##End while eps > .00001)


#-------Compute ML estimate of standard error-(Hamdan, 1970)----#
#  Hamdan, M. A.  (1970).  The Equivalence of Tetrachoric and Maximum 
#  Likelihood Estimates of r in 2 x 2 Tables.  Biometrika, 57(1), 212-215.  
 if( stderror ){
   se[iter.row,iter.col] <- 1/(N * bvn(c(h.x,k.y))) *  (1/cella + 1/cellb + 1/cellc + 1/celld)^-.5
 } ## END 0F if(stderror)
#-------------------------------------------------------#

   r[iter.row,iter.col]<-rhat 
   }## END iter over columns
  } ## END iter over rows

## fill in upper triangle of r 
 r <- r + t(r)
 diag(r) <- 1

## Smooth tetrachoric R matrix
 if(Smooth){
     ULU <- eigen(r)
     U <- ULU$vectors
     L <- ULU$values

     
     
     ## maybe use smoothAPA here
     if( min(L) <= 0 ){  # renorm to make matrix positive definite     
         L[L<=0] <- .0001
         Ltot <- sum(L)
         L <- nitems*L/Ltot
         Lsqrt <- diag(sqrt(L))
         Fload <- tcrossprod(U, Lsqrt)       #U%*%Lsqrt
         r <- tcrossprod(Fload, Fload)     #Fload %*% t(Fload)
         Dmat <- diag(1/sqrt(diag(r)))
         r <- Dmat %*% r %*% Dmat
     }
 }##End if Smooth   



## Return results
 if( stderror == FALSE ) {
      list(r = r, 
           se = NA,
           convergence = converged,
           Warnings = warning.msg)
    }

 else if ( stderror == TRUE ){
       se <- se + t(se)
       diag(se) <- NA
       ## If a standard error cannot be computed it will be reported as 
       ## NA.
       se[seFlagMatrix==99] <-NA
       list(r = r, 
            se = se,
            convergence = converged,
            Warnings = warning.msg)
  } 

} ## End of tetcor function



######### CHECK #########
# ## generate bivariate normal data
# library(MASS)
# library(microbenchmark)
# rho <- .85
# xy<- mvrnorm(50000, mu = c(0,0), Sigma = matrix(c(1,rho,rho,1),ncol=2))
# 
# ## dichotomize
# p1 <- .59
# p2 <- .1
# xy[,1] <- (xy[,1] > p1) 
# xy[,2] <- (xy[,2] > p2)
# 
# 
# 
# tetcor(X=xy[,1], y=xy[,2], max.iter=5000,
#        stderror = TRUE, PRINT = FALSE)
# 
# 
# microbenchmark(times = 200,
# tetcor(X=xy[,1], y=xy[,2], max.iter=5000,
#        stderror = TRUE, PRINT = FALSE)
# )
# 


