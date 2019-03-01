## FMP 
# programmed by Niels Waller 
# January 27, 2016
# bugs in T2 T3 fixed on February 8, 2016



#########################################################################


#' Estimate the coefficients of a filtered monotonic polynomial IRT model
#' 
#' Estimate the coefficients of a filtered monotonic polynomial IRT model.
#' 
#' As described by Liang and Browne (2015), the filtered polynomial model (FMP)
#' is a quasi-parametric IRT model in which the IRF is a composition of a
#' logistic function and a polynomial function, \eqn{m(\theta)}, of degree 2k +
#' 1.  When k = 0, \eqn{m(\theta) = b_0 + b_1 \theta} (the slope intercept form
#' of the 2PL). When k = 1, 2k + 1 equals 3 resulting in \eqn{m(\theta) = b_0 +
#' b_1 \theta + b_2 \theta^2 + b_3 \theta^3}. Acceptable values of k = 0,1,2,3.
#' According to Liang and Browne, the "FMP IRF may be used to approximate any
#' IRF with a continuous derivative arbitrarily closely by increasing the
#' number of parameters in the monotonic polynomial" (2015, p. 2) The FMP model
#' assumes that the IRF is monotonically increasing, bounded by 0 and 1, and
#' everywhere differentiable with respect to theta (the latent trait).
#' 
#' @param data N(subjects)-by-p(items) matrix of 0/1 item response data.
#' @param thetaInit Initial theta (\eqn{\theta}) surrogates (e.g., calculated
#' by \link{svdNorm}).
#' @param item Item number for coefficient estimation.
#' @param startvals Start values for function minimization. Start values are in
#' the gamma metric (see Liang & Browne, 2015)
#' @param k Order of monotonic polynomial = 2k+1 (see Liang & Browne, 2015). k
#' can equal 0, 1, 2, or 3.
#' @param eps Step size for gradient approximation, default = 1e-6. If a
#' convergence failure occurs during function optimization reducing the value
#' of eps will often produce a converged solution.
#' @return \item{b}{Vector of polynomial coefficients.} \item{gamma}{Polynomial
#' coefficients in gamma metric (see Liang & Browne, 2015).}
#' \item{FHAT}{Function value at convergence.} \item{counts}{Number of function
#' evaluations during minimization (see optim documentation for further
#' details).} \item{AIC}{Pseudo scaled Akaike Information Criterion (AIC).
#' Candidate models that produce the smallest AIC suggest the optimal number of
#' parameters given the sample size. Scaling is accomplished by dividing the
#' non-scaled AIC by sample size.} \item{BIC}{Pseudo scaled Bayesian
#' Information Criterion (BIC). Candidate models that produce the smallest BIC
#' suggest the optimal number of parameters given the sample size. Scaling is
#' accomplished by dividing the non-scaled BIC by sample size.}
#' \item{convergence}{Convergence = 0 indicates that the optimization algorithm
#' converged; convergence=1 indicates that the optimization failed to
#' converge.}
#' @author Niels Waller
#' @references Liang, L. & Browne, M. W. (2015). A quasi-parametric method for
#' fitting flexible item response functions. \emph{Journal of Educational and
#' Behavioral Statistics, 40}, 5--34.
#' @keywords statistics
#' @export
#' @examples
#' 
#' 
#' \dontrun{
#' ## In this example we will generate 2000 item response vectors 
#' ## for a k = 1 order filtered polynomial model and then recover 
#' ## the estimated item parameters with the FMP function.  
#' 
#' k <- 1  # order of polynomial
#' 
#' NSubjects <- 2000
#' 
#' 
#' ## generate a sample of 2000 item response vectors 
#' ## for a k = 1 FMP model using the following
#' ## coefficients
#' b <- matrix(c(
#'    #b0     b1      b2     b3   b4  b5  b6  b7  k
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
#' ex1.data<-genFMPData(NSubj = NSubjects, bParams = b, seed = 345)$data
#' 
#' ## number of items in the data matrix
#' NItems <- ncol(ex1.data)
#' 
#' # compute (initial) surrogate theta values from 
#' # the normed left singular vector of the centered 
#' # data matrix
#' thetaInit <- svdNorm(ex1.data)
#' 
#' 
#' ## earlier we defined k = 1
#'   if(k == 0) {
#'             startVals <- c(1.5, 1.5)
#'             bmat <- matrix(0, NItems, 6)
#'             colnames(bmat) <- c(paste("b", 0:1, sep = ""),"FHAT", "AIC", "BIC", "convergence") 
#'   }
#'   if(k == 1) {
#'            startVals <- c(1.5, 1.5, .10, .10)
#'            bmat <- matrix(0, NItems, 8)
#'            colnames(bmat) <- c(paste("b", 0:3, sep = ""),"FHAT", "AIC", "BIC", "convergence") 
#'   }
#'   if(k == 2) {
#'            startVals <- c(1.5, 1.5, .10, .10, .10, .10)
#'            bmat <- matrix(0, NItems, 10)
#'            colnames(bmat) <- c(paste("b", 0:5, sep = ""),"FHAT", "AIC", "BIC", "convergence") 
#'   }
#'   if(k == 3) {
#'            startVals <- c(1.5, 1.5, .10, .10, .10, .10, .10, .10)
#'            bmat <- matrix(0, NItems, 12)
#'            colnames(bmat) <- c(paste("b", 0:7, sep = ""),"FHAT", "AIC", "BIC", "convergence") 
#'   }         
#'   
#' # estimate item parameters and fit statistics  
#'   for(i in 1:NItems){
#'     out <- FMP(data = ex1.data, thetaInit, item = i, startvals = startVals, k = k)
#'     Nb <- length(out$b)
#'     bmat[i,1:Nb] <- out$b
#'     bmat[i,Nb+1] <- out$FHAT
#'     bmat[i,Nb+2] <- out$AIC
#'     bmat[i,Nb+3] <- out$BIC
#'     bmat[i,Nb+4] <- out$convergence
#'   }
#' 
#' # print output 
#' print(bmat)
#' }
#' 
FMP<-function(data, thetaInit, item, startvals, k=0, eps=1e-6){
  
	if(k > 3) stop("\nFatal Error: k must equal 0, 1, 2, or 3.\n")
  
    T1 <- matrix(0,3,1)
    T2 <- matrix(0,5,3); 
    T2[1,1]<-T2[2,2]<-T2[3,3]<-1
    T3 <- matrix(0,7,5)
    T3[1,1]<-T3[2,2]<-T3[3,3]<-T3[4,4]<-T3[5,5]<-1
   
  
    # Prob of a keyed response 2PL
    P <- function(m){
        1/(1+exp(-m))
    }  
  
    #------------------------------------------------------------------##
    # gammaTob:  generate final b coefficients from final gamma estimates
    gammaTob<-function(gammaParam, k=1){
    
           if(k>3) stop("\n\n*** FATAL ERROR:k must be <=3 ***\n\n")
    
           if(k==0){
             Ksi = gammaParam[1]
             lambda = gammaParam[2]
             b0 = Ksi
             b1 = lambda
             return(c(b0,b1))
           } 
    
           if(k==1){
             Ksi = gammaParam[1]
             lambda = gammaParam[2]
             alpha1<-gammaParam[3]
             beta1 <-gammaParam[4]
             phi1<-alpha1^2 + beta1
             T1 <-as.vector(c(1, -2*alpha1, phi1))
             asup0 = lambda
             asup1 = T1 * asup0
             b0 = Ksi
             b1 = asup1[1]/1
             b2 = asup1[2]/2
             b3 = asup1[3]/3
            return( c(b0,b1,b2,b3))
          }  
    
          if(k==2){
            Ksi = gammaParam[1]
            lambda = gammaParam[2]
            alpha1 <- gammaParam[3]
            beta1  <- gammaParam[4]
            alpha2 <- gammaParam[5]
            beta2  <- gammaParam[6]
            phi1   <- alpha1^2 + beta1
            phi2   <- alpha2^2 + beta2
            T1 <-as.vector(c(1, -2*alpha1, phi1))
            T2[2,1]<-T2[3,2]<-T2[4,3] <- -2*alpha2
            T2[3,1] <- T2[4,2] <- T2[5,3]<- phi2
            # see (A6) Liang & Browne
            asup1 = (T2 %*% T1) * lambda
            b0 = Ksi
            b1 = asup1[1]/1
            b2 = asup1[2]/2
            b3 = asup1[3]/3
            b4 = asup1[4]/3
            b5 = asup1[5]/3
            return( c(b0,b1,b2,b3,b4,b5))
          }  
    
          if(k==3){
            Ksi = gammaParam[1]
            lambda = gammaParam[2]
      
            alpha1 <- gammaParam[3]
            beta1  <- gammaParam[4]
            alpha2 <- gammaParam[5]
            beta2  <- gammaParam[6]
            alpha3 <- gammaParam[7]
            beta3  <- gammaParam[8]
      
            phi1 <- alpha1^2 + beta1
            phi2 <- alpha2^2 + beta2
            phi3 <- alpha3^2 + beta3
      
            T1 <-as.vector(c(1, -2*alpha1, phi1))
            T2[2,1]<-T2[3,2]<-T2[4,3] <- -2*alpha2
            T2[3,1] <- T2[4,2] <- T2[5,3]<- phi2
            T3[2,1]<-T3[3,2]<-T3[4,3]<-T3[5,4]<-T3[6,5]<- -2*alpha3
            T3[3,1]<-T3[4,2]<-T3[5,3]<-T3[6,4]<-T3[7,5]<- phi3
      
            # see (A6) Liang & Browne
            asup3 = (T3 %*% T2 %*% T1) * lambda
      
            b0 = Ksi
            b1 = asup3[1]/1
            b2 = asup3[2]/2
            b3 = asup3[3]/3
            b4 = asup3[4]/4
            b5 = asup3[5]/5
            b6 = asup3[6]/6
            b7 = asup3[7]/7
            return( c(b0,b1,b2,b3,b4,b5,b6,b7))
         }  
  } # END gammaTob -----------------------------------------------##
  
  
  
  
  
  ###################################################################
  ## fncFMP
    fncFMP<- function(gammaParam, k=1, data, item, thetaInit){
    
          if(k>3) stop("\n\n*** FATAL ERROR:k must be <=3 ***\n\n")
    
          if(k==0){
                Ksi = gammaParam[1]
                lambda = gammaParam[2]
         
                b0 = Ksi
                b1 = lambda
      
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
                Ksi = gammaParam[1]
                lambda = gammaParam[2]
      
                alpha1<-gammaParam[3]
                beta1 <-gammaParam[4]
      
                phi1<-alpha1^2 + beta1
      
               T1 <-as.vector(c(1, -2*alpha1, phi1))
      
               asup0 = lambda
               asup1 = T1 * asup0
               b0 = Ksi
               b1 = asup1[1]/1
               b2 = asup1[2]/2
               b3 = asup1[3]/3
      
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
                Ksi = gammaParam[1]
                lambda = gammaParam[2]
      
                alpha1 <- gammaParam[3]
                beta1  <- gammaParam[4]
                alpha2 <- gammaParam[5]
                beta2  <- gammaParam[6]
      
                phi1 <- alpha1^2 + beta1
                phi2 <- alpha2^2 + beta2
      
                T1 <-as.vector(c(1, -2*alpha1, phi1))
                T2[2,1]<-T2[3,2]<-T2[4,3] <- -2*alpha2
                T2[3,1] <- T2[4,2] <- T2[5,3]<- phi2
      
                # see (A6) Liang & Browne
                asup2 = (T2 %*% T1) * lambda
      
                b0 = Ksi
                b1 = asup2[1]/1
                b2 = asup2[2]/2
                b3 = asup2[3]/3
                b4 = asup2[4]/4
                b5 = asup2[5]/5
      
                y = data[,item]
      
                # create logit polynomial 
                m = b0 + b1*thetaInit + b2*thetaInit^2 + b3*thetaInit^3 +    b4*thetaInit^4 + b5*thetaInit^5
      
                ##  avoid log(0)
                Pm <- P(m);
      
                Pm[Pm=="NaN"]<-.5
                Pm[Pm==1]   <-1-1e-12
                Pm[Pm<1e-12]<-1e-12  
      
                # this is FHAT
                return(-mean( y*log(Pm) + (1 - y) * log(1 - Pm)  ))
          }#END k=2 
    
    
          if(k==3){
                Ksi = gammaParam[1]
                lambda = gammaParam[2]
      
                alpha1 <- gammaParam[3]
                beta1  <- gammaParam[4]
                alpha2 <- gammaParam[5]
                beta2  <- gammaParam[6]
                alpha3 <- gammaParam[7]
                beta3  <- gammaParam[8]
      
                phi1 <- alpha1^2 + beta1
                phi2 <- alpha2^2 + beta2
                phi3 <- alpha3^2 + beta3
      
                T1 <-as.vector(c(1, -2*alpha1, phi1))
                T2[2,1]<-T2[3,2]<-T2[4,3] <- -2*alpha2
                T2[3,1] <- T2[4,2] <- T2[5,3]<- phi2
                T3[2,1]<-T3[3,2]<-T3[4,3]<-T3[5,4]<-T3[6,5]<- -2*alpha3
                T3[3,1]<-T3[4,2]<-T3[5,3]<-T3[6,4]<-T3[7,5]<- phi3
                # see (A6) Liang & Browne
                asup3 = (T3 %*% T2 %*% T1) * lambda
      
                b0 = Ksi
                b1 = asup3[1]/1
                b2 = asup3[2]/2
                b3 = asup3[3]/3
                b4 = asup3[4]/4
                b5 = asup3[5]/5
                b6 = asup3[6]/6
                b7 = asup3[7]/7
      
      
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
    
  } #END fncFMP ##################################################

  
  
  
  
  if(k>3) stop("\n\n*** FATAL ERROR:k must be <= 3 ***\n\n")
  
    if(k==0) Constraints = c(-Inf,0)
    if(k==1) Constraints = c(-Inf,0,-Inf,0)
    if(k==2) Constraints = c(-Inf,0,-Inf,0,-Inf,0)
    if(k==3) Constraints = c(-Inf,0,-Inf,0,-Inf,0,-Inf,0)
  
    Nparam<-length(startvals)
    
    out <- optim(par=startvals, 
                 fn=fncFMP, 
                 method="L-BFGS-B",
                 lower= Constraints,
                 k=k,
                 data=data,
                 item=item,
                 thetaInit=thetaInit,
                 control=list(ndeps=rep(eps,Nparam),
                              maxit=2000))
    
  ## compute scaled (pseudo) AIC and BIC
     NSubj <- nrow(data)
     q <- 2 * k + 2
     AIC <- 2 * out$value + (2/NSubj) * q
     BIC <- 2 * out$value + (log(NSubj)/NSubj) * q
  
     
   
     list( b = gammaTob(out$par,k),
         gamma=out$par, 
         FHAT=out$value,
         counts=out$counts,
         AIC = AIC,
         BIC = BIC,
         convergence = out$convergence)
  
} #END FMP




###################################################################
#                 END OF FUNCTION DEFINITIONS
###################################################################



