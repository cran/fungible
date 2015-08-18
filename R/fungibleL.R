#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# fungibleL:  A program to compute Fungible Logistic Regression Coefficients
# December 13, 2013
# Jeff Jones and Niels Waller
# args:
#  X          : an N by Nvar matrix of predictor scores without the leading 
#               column of ones
#  y          : and N by 1 vector of dichotomous criterion scores
#  Nsets      : the desired number of fungible coefficient vectors
#  method     : LLM: Log-Likelihood method. EM: Ellipsoid Method
#               default: method="LLM"
#  RsqDelta   : the desired decrement in the pseudo-R-squared
#  rLaLb      : the desired correlation between the logits
#  s          : scale factor for random deviates
# 
# output:
#  model      : the glm model object
#  call       : the function call
#  ftable     : a data frame with the mle estimates and the minimum and maximum fungible
#               coefficients
#  lnLML      : the maximum likelihood log likelihood value
#  lnLf       : the decremented, fungible log likelihood value
#  pseudoRsq  : the pseudo R-squared
#  fungibleRsq: the fungible pseudo R-squared
#  fungibleb  : the Nsets by Nvar + 1 matrix of fungible coefficients
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  

fungibleL <- function(X, y, Nsets = 1000, method = "LLM", RsqDelta = 0.01, rLaLb = NULL, s = .3) {
  
  #library(MASS)
  
  Xtmp <- X
  X    <- cbind(1, X)
  npar <- ncol(X); N <- nrow(X)

  fungibleb <- matrix(0, Nsets, npar)
  
  mod <- glm(y ~ Xtmp, family="binomial")
  lnLML <- as.numeric(logLik(mod))
  # maximum likelihood coefficient estimates
  bmle <- coef(mod)
  
  cl = match.call()
  
  if(method=="LLM"){
  
  # Log-Likelihood Function
    ll <- function(b, y, X, ct.pt) {
      sum(y*log((1+exp(-X%*%b))^-1)) +  sum( (1-y) * log( (1 - (1+exp(-X%*%b))^-1))) - ct.pt
    }     
  
  # First Derivative of Log-Likelihood function
    dll1 <- function(b, y, X){
      Z <- X %*% b
      P <- 1/(1 + exp(-Z))
      apply(  as.numeric(y*(1-P))*X - as.numeric((1-y)*P)*X, 2, sum)
    }
  
  # log likelihood ratio pseudo R-squared
    RsqFUN <- 1 - mod$dev/mod$null
  
    if(RsqFUN - RsqDelta < 0) 
      stop("Model Decrement out of Bounds. Select a smaller RsqDelta.\n")

  # height of log-likelihood for fungible weights
    lnLf <- (-mod$null/2)*(1-(RsqFUN-RsqDelta))  
      
  # a matrix of random start values  
    sVcovMod<-s*vcov(mod)
    rstarts <- matrix(MASS::mvrnorm(n=Nsets,mu=bmle,Sigma=sVcovMod), Nsets, npar, byrow=T)

  # create a progress bar 
    pb <- txtProgressBar(min=1, max=Nsets, style=3)
    cat("\n   Computing fungible regression coefficients . . .\n")
      
    for(i in 1:Nsets){
      
      tol <- 99
      bi <- rstarts[i,]
      
      while(tol > 1e-8) {
      
        g <- dll1(bi, y, X)                  
        gInv <- solve(t(g)%*%g)*g
            
      # Newton Root Finding       
        biPlus1 <- bi - gInv * ll(bi, y, X, ct.pt=lnLf)
                  
        if(is.na(ll(biPlus1, y, X, ct.pt=lnLf))) {         
          
          bi <- MASS::mvrnorm(n=1,mu=bmle,Sigma=sVcovMod)
        
        } else { 
          
          bi <- biPlus1
          tol <- abs(ll(biPlus1, y, X, ct.pt=lnLf))
        
        }
      
      }# END while(tol > 1e-8) 
    
      fungibleb[i,] <- biPlus1 
      setTxtProgressBar(pb, i)
  
    } #END for(i in 1:Nsets)

    close(pb)  
    
  # summarize results  
    fmin <- apply(fungibleb,2,min)
    fmax <- apply(fungibleb,2,max)
    ftable = data.frame("mle estimates" = round(coef(mod),3),
                        "fungible.min"=round(fmin,3),
                        "fungible.max"=round(fmax,3))

    rownames(ftable) <- stringr::str_replace(rownames(ftable), "Xtmp", "")

    list(model=mod, call = cl, ftable = ftable, fungiblea = fungibleb, 
         lnLML=lnLML,lnLf=lnLf, "pseudoRsq" = RsqFUN,
         rLaLb=NULL, "fungibleRsq"=RsqFUN-RsqDelta)  
    
  } else { # Begin ellipsoid method
   
    X <- X[,-1]
    
    # norms the columns of a matrix  
    norm <- function(x) {			
      x <- as.matrix(x)
      xx <- t(x)%*%x
      x%*%diag(diag(xx)^-.5,nrow=nrow(xx),ncol=ncol(xx))
    }
    
    log.odds <- log(mod$fitted/(1-mod$fitted))
    ymn <- mean(log.odds); Xmn <- apply(X, 2, mean)
    
    Sxx <- cov(X); sxy <- cov(X,log.odds)
            
    rLaLb <- rLaLb 
    
    p <- ncol(Sxx)
      
    UDU <- eigen(Sxx)
      
    U <- UDU$vector
    D <- diag(UDU$val)
    Dhalf <- sqrt(D)
    Dinvhalf <- diag(1/diag(Dhalf))
          
    g <- sqrt(D)%*%t(U)%*%bmle[-1]/as.numeric(sqrt(t(bmle[-1])%*%Sxx%*%bmle[-1]))
              
  # create a progress bar 
    pb <- txtProgressBar(min=1, max=Nsets, style=3)
    cat("\n   Computing fungible regression coefficients . . .\n")
    
  # Generate fungible weights 
      
    for(i in 1:Nsets) {
        
      tmp <- matrix(c(g,rnorm(p*(p-1))),p,p)
      O <- qr.Q(qr(tmp))[,-1]    
        
      v <- norm(rnorm((p-1)))
        
      hn <- rLaLb*g + sqrt(1-rLaLb^2)*O%*%v
        
      atil <- U%*%Dinvhalf%*%hn
        
      s <- as.numeric(t(sxy)%*%atil)
      s.sign <- sign(s)
        
      fungibleb[i,-1] <- s.sign*s*atil	        
      fungibleb[i,1] <- ymn - Xmn%*%fungibleb[i,-1]         
    
      setTxtProgressBar(pb, i)
      
    } # end ellipsoid method
    
    close(pb)
    
    # summarize results  
    fmin <- apply(fungibleb,2,min)
    fmax <- apply(fungibleb,2,max)
    ftable = data.frame("mle estimates" = round(coef(mod),3),
                        "fungible.min"=round(fmin,3),
                        "fungible.max"=round(fmax,3))

    rownames(ftable) <- stringr::str_replace(rownames(ftable), "Xtmp", "")

    list(model=mod, call = cl, ftable = ftable, fungiblea = fungibleb, 
         lnLML=NULL,lnLf=NULL, "pseudoRsq" = NULL,
         rLaLb=rLaLb, "fungibleRsq"=NULL)
    
  } # end else statement 
        
}# END of function fungibleL



