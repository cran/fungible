## FUP in R
# programmed by Niels Waller 
# October 14, 2014

###################################################################


#' Estimate the coefficients of a filtered unconstrained polynomial IRT model
#' 
#' Estimate the coefficients of a filtered unconstrained polynomial IRT model.
#' 
#' 
#' @param data N(subjects)-by-p(items) matrix of 0/1 item response data.
#' @param thetaInit Initial theta surrogates (e.g., calculated by
#' \link{svdNorm}).
#' @param item item number for coefficient estimation.
#' @param startvals start values for function minimization.
#' @param k order of monotonic polynomial = 2k+1 (see Liang & Browne, 2015).
#' @return \item{b}{Vector of polynomial coefficients.} \item{FHAT}{Function
#' value at convergence.} \item{counts}{Number of function evaluations during
#' minimization (see optim documentation for further details).}
#' \item{AIC}{Pseudo scaled Akaike Information Criterion (AIC). Candidate
#' models that produce the smallest AIC suggest the optimal number of
#' parameters given the sample size. Scaling is accomplished by dividing the
#' non-scaled AIC by sample size.} \item{BIC}{Pseudo scaled Bayesian
#' Information Criterion (BIC). Candidate models that produce the smallest BIC
#' suggest the optimal number of parameters given the sample size. Scaling is
#' accomplished by dividing the non-scaled BIC by sample size.}
#' \item{convergence}{Convergence = 0 indicates that the optimization algorithm
#' converged; convergence=1 indicates that the optimization failed to converge.
#' 
#' .}
#' @author Niels Waller
#' @references Liang, L. & Browne, M. W. (2015). A quasi-parametric method for
#' fitting flexible item response functions. \emph{Journal of Educational and
#' Behavioral Statistics, 40}, 5--34.
#' @keywords statistics
#' @export
#' @examples
#' 
#' \dontrun{
#' NSubjects <- 2000
#' 
#' 
#' ## generate sample k=1 FMP data
#' b <- matrix(c(
#'     #b0    b1     b2    b3      b4   b5 b6 b7  k
#'   1.675, 1.974, -0.068, 0.053,  0,  0,  0,  0, 1,
#'   1.550, 1.805, -0.230, 0.032,  0,  0,  0,  0, 1,
#'   1.282, 1.063, -0.103, 0.003,  0,  0,  0,  0, 1,
#'   0.704, 1.376, -0.107, 0.040,  0,  0,  0,  0, 1,
#'   1.417, 1.413,  0.021, 0.000,  0,  0,  0,  0, 1,
#'  -0.008, 1.349, -0.195, 0.144,  0,  0,  0,  0, 1,
#'   0.512, 1.538, -0.089, 0.082,  0,  0,  0,  0, 1,
#'   0.122, 0.601, -0.082, 0.119,  0,  0,  0,  0, 1,
#'   1.801, 1.211,  0.015, 0.000,  0,  0,  0,  0, 1,
#'  -0.207, 1.191,  0.066, 0.033,  0,  0,  0,  0, 1,
#'  -0.215, 1.291, -0.087, 0.029,  0,  0,  0,  0, 1,
#'   0.259, 0.875,  0.177, 0.072,  0,  0,  0,  0, 1,
#'  -0.423, 0.942,  0.064, 0.094,  0,  0,  0,  0, 1,
#'   0.113, 0.795,  0.124, 0.110,  0,  0,  0,  0, 1,
#'   1.030, 1.525,  0.200, 0.076,  0,  0,  0,  0, 1,
#'   0.140, 1.209,  0.082, 0.148,  0,  0,  0,  0, 1,
#'   0.429, 1.480, -0.008, 0.061,  0,  0,  0,  0, 1,
#'   0.089, 0.785, -0.065, 0.018,  0,  0,  0,  0, 1,
#'  -0.516, 1.013,  0.016, 0.023,  0,  0,  0,  0, 1,
#'   0.143, 1.315, -0.011, 0.136,  0,  0,  0,  0, 1,
#'   0.347, 0.733, -0.121, 0.041,  0,  0,  0,  0, 1,
#'  -0.074, 0.869,  0.013, 0.026,  0,  0,  0,  0, 1,
#'   0.630, 1.484, -0.001, 0.000,  0,  0,  0,  0, 1), 
#'   nrow=23, ncol=9, byrow=TRUE)  
#'  
#' # generate data using the above item parameters 
#' ex1.data<-genFMPData(NSubj = NSubjects, bParams = b, seed = 345)$data
#' 
#' NItems <- ncol(ex1.data)
#' 
#' # compute (initial) surrogate theta values from 
#' # the normed left singular vector of the centered 
#' # data matrix
#' thetaInit <- svdNorm(ex1.data)
#' 
#' # Choose model
#' k <- 1  # order of polynomial = 2k+1
#' 
#' # Initialize matrices to hold output
#' if(k == 0) {
#'   startVals <- c(1.5, 1.5)
#'   bmat <- matrix(0,NItems,6)
#'   colnames(bmat) <- c(paste("b", 0:1, sep = ""),"FHAT", "AIC", "BIC", "convergence") 
#' }
#' 
#' if(k == 1) {
#'   startVals <- c(1.5, 1.5, .10, .10)
#'   bmat <- matrix(0,NItems,8)
#'   colnames(bmat) <- c(paste("b", 0:3, sep = ""),"FHAT", "AIC", "BIC", "convergence") 
#' }
#' 
#' if(k == 2) {
#'   startVals <- c(1.5, 1.5, .10, .10, .10, .10)
#'   bmat <- matrix(0,NItems,10)
#'   colnames(bmat) <- c(paste("b", 0:5, sep = ""),"FHAT", "AIC", "BIC", "convergence") 
#' }
#' 
#' if(k == 3) {
#'   startVals <- c(1.5, 1.5, .10, .10, .10, .10, .10, .10)
#'   bmat <- matrix(0,NItems,12)
#'   colnames(bmat) <- c(paste("b", 0:7, sep = ""),"FHAT", "AIC", "BIC", "convergence") 
#' }   
#' 
#' 
#' # estimate item parameters and fit statistics
#' for(i in 1:NItems){
#'   out<-FUP(data = ex1.data,thetaInit = thetaInit, item = i, startvals = startVals, k = k)
#'   Nb <- length(out$b)
#'   bmat[i,1:Nb] <- out$b
#'   bmat[i,Nb+1] <- out$FHAT
#'   bmat[i,Nb+2] <- out$AIC
#'   bmat[i,Nb+3] <- out$BIC
#'   bmat[i,Nb+4] <- out$convergence
#' }
#' 
#' # print results
#' print(bmat)
#' }
#' 
FUP<-function(data, thetaInit, item, startvals, k = 0){
  
  
  # Prob of a keyed response 2PL ------------------------------##
  P <- function(m){
    1/(1+exp(-m))
  }
  
  ##-----------------------------------------------------------##
  fncFUP<- function(bParam, k=1, data, item, thetaInit){
    
       if(k>3) stop("\n\n*** FATAL ERROR:k must be <=3 ***\n\n")
    
       if(k==0){
      
         b0 = bParam[1]
         b1 = bParam[2]
      
         y = data[,item]
      
         m = b0 + b1*thetaInit 
      
         ##  avoid log(0)
         Pm <- P(m);
         Pm[Pm==1]   <-1-1e-12
         Pm[Pm<1e-12]<-1e-12
      
         # this is FHAT
         return(-mean( y*log(Pm) + (1 - y) * log(1 - Pm)  ))
      }#END k=0 
    
      if(k==1){
      
        b0 = bParam[1]
        b1 = bParam[2]
        b2 = bParam[3]
        b3 = bParam[4]
      
        y = data[,item]
      
        m = b0 + b1*thetaInit + b2*thetaInit^2 + b3*thetaInit^3
      
        ##  avoid log(0)
        Pm <- P(m);
      
        Pm[Pm=="NaN"]<-.5
        Pm[Pm==1]   <-1-1e-12
        Pm[Pm<1e-12]<-1e-12  
      
        # this is FHAT
        return(-mean( y*log(Pm) + (1 - y) * log(1 - Pm)  ))
      }#END k=1 
    
      if(k==2){
      
        b0 = bParam[1]
        b1 = bParam[2]
        b2 = bParam[3]
        b3 = bParam[4]
        b4 = bParam[5]
        b5 = bParam[6]
      
        y = data[,item]
      
        # create logit polynomial 
        m = b0 + b1*thetaInit + b2*thetaInit^2 + b3*thetaInit^3 + b4*thetaInit^4 + b5*thetaInit^5
      
        ##  avoid log(0)
        Pm <- P(m);
      
        Pm[Pm=="NaN"]<-.5
        Pm[Pm==1]   <-1-1e-12
        Pm[Pm<1e-12]<-1e-12  
      
        # this is FHAT
        return(-mean( y*log(Pm) + (1 - y) * log(1 - Pm)  ))
      }#END k=2 
    
    
      if(k==3){
      
        b0 = bParam[1]
        b1 = bParam[2]
        b2 = bParam[3]
        b3 = bParam[4]
        b4 = bParam[5]
        b5 = bParam[6]
        b6 = bParam[7]
        b7 = bParam[8]
      
        y = data[,item]
      
        # create logit polynomial 
        m = b0 + b1*thetaInit + b2*thetaInit^2 + b3*thetaInit^3 + 
            b4*thetaInit^4 + b5*thetaInit^5    +
            b6*thetaInit^6 + b7*thetaInit^7
      
        ##  avoid log(0)
        Pm <- P(m);
      
        Pm[Pm=="NaN"]<-.5
        Pm[Pm==1]   <-1-1e-12
        Pm[Pm<1e-12]<-1e-12  
      
        # this is FHAT
        return(-mean( y*log(Pm) + (1 - y) * log(1 - Pm)  ))
      }#END k=3 
    
  } #END fncFUP --------------------------------------------##
  
 
  
    out <- optim(par=startvals, 
                 fn=fncFUP, 
                 method="BFGS",
                 k=k,
                 data=data,
                 item=item,
                 thetaInit=thetaInit,
                 control=list(factr=1e-12,
                              maxit=1000))
    
  ## compute scaled (pseudo) AIC and BIC
     NSubj <- nrow(data)
     q <- 2 * k + 2
     AIC <- 2 * out$value + (2/NSubj) * q
     BIC <- 2 * out$value + (log(NSubj)/NSubj) * q
  
 
   list( b = out$par,
        FHAT=out$value,
        counts=out$counts,
        AIC = AIC,
        BIC = BIC,
        convergence = out$convergence)
  
} #END FMP




###################################################################
#                 END OF FUNCTION DEFINITIONS
###################################################################



