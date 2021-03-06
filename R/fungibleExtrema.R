 # December 31, 2020 updated argument list
# Version: October 23, 2012
# August 6, 2018 (updated recycling)
#                    R Code for fungibleExtrema
#
##  R function: fungibleExtrema
##
## Locate extrema of Fungible weights (i.e., find max (min) cos(b,a)
##  
##
##  Input Variables
##  R.X   p x p Predictor variable correlation matrix.
##  rxy   p x 1 Vector of predictor-criterion correlations.
##  r.yhata.yhatb  = correlation between alternate-weight (yhata) and 
##                   and  least squares (yhatb) composites.
##  Nstarts  Maximum number of (max) minimizations from random starting configurations.
##  MaxMin   Max = maximize cos(a,b); Min = minimize cos(a,b)
##
##  Output Variables
##  cos.ab  cosine between OLS and alternate weights
##  a     extrema of fungible weights
##  k     k weights
##  z     z weights
##  b     OLS weights
##  u     u weights
##  r.yhata.yhatb  correlation between yhat and ytilde
##  r.y.yhat       correlation between y and yhat
##  gradient   gradient at solution

        



#' Locate Extrema of Fungible Regression Weights
#' 
#' Locate extrema of fungible regression weights.
#' 
#' 
#' @param R.X p x p Predictor variable correlation matrix.
#' @param rxy p x 1 Vector of predictor-criterion correlations.
#' @param r.yhata.yhatb Correlation between least squares (yhatb) and
#' alternate-weight (yhata) composites.
#' @param Nstarts Maximum number of (max) minimizations from random starting
#' configurations.
#' @param MaxMin Character: "Max" = maximize cos(a,b); "Min" = minimize
#' cos(a,b).
#' @param Seed Starting seed for the random number generator. If Seed = NULL 
#' then the program will sample a random integer in the (0, 100,000) interval.
#' Default (Seed = NULL).  
#' @param maxGrad The optimization routine will end when the maximimum of
#' the (absolute value of the ) function gradient falls below the value specified in 
#' maxGrad. Default (maxGrad = 1E-05). 
#' @param PrintLevel (integer). If PrintLevel = 1 then the program will print 
#' additional output during function convergence. Default (PrintLevel = 1).
#' @return \item{cos.ab}{cosine between OLS and alternate weights.}
#' \item{a}{extrema of fungible weights.} \item{k}{k weights.} \item{z}{z
#' weights: a normalized random vector.} \item{b}{OLS weights.} \item{u}{p x 1
#' vector of u weights.} \item{r.yhata.yhatb}{Correlation between yhata and
#' yhatb.} \item{r.y.yhatb}{Correlation between y and yhatb.}
#' \item{gradient}{Gradient of converged solution.}
#' 
#' @author Niels Waller and Jeff Jones
#' 
#' @references Koopman, R. F.  (1988).  On the sensitivity of a composite to
#' its weights.  \emph{Psychometrika, 53(4)}, 547--552.
#' 
#' Waller, N. & Jones, J. (2009). Locating the extrema of fungible regression
#' weights in multiple regression. \emph{Psychometrika, 74}, 589--602.
#' @keywords fungible
#' @export
#' @examples
#' \dontrun{  
#' ## Example 
#' ## This is Koopman's Table 2 Example
#' 
#' 
#' R.X <- matrix(c(1.00,  .69,  .49,  .39,
#'                 .69, 1.00,  .38,  .19,
#'                 .49,  .38, 1.00,  .27,
#'                 .39,  .19,  .27, 1.00),4,4)
#' 
#' 
#' b <- c(.39, .22, .02, .43)
#' rxy <- R.X %*% b
#' 
#' OLSRSQ <- t(b) %*% R.X %*% b
#' 
#' theta <- .02
#' r.yhata.yhatb <- sqrt( 1 - (theta)/OLSRSQ)
#' 
#' Converged = FALSE
#' SEED = 1234
#' MaxTries = 100 
#' iter = 1
#' 
#' while( iter <= MaxTries){
#'    SEED <- SEED + 1
#'   
#'    cat("\nCurrent Seed = ", SEED, "\n")
#'    output <- fungibleExtrema(R.X, rxy, 
#'                              r.yhata.yhatb, 
#'                              Nstarts = 5,
#'                              MaxMin = "Min", 
#'                              Seed = SEED,
#'                              maxGrad = 1E-05,
#'                              PrintLevel = 1)
#'   
#'    Converged <- output$converged
#'    if(Converged) break
#'    iter = iter + 1
#' }  
#' 
#' print( output )
#' 
#' ## Scale to replicate Koopman
#' a <- output$a
#' a.old <- a
#' aRa <- t(a) %*% R.X %*% a
#' 
#' ## Scale a such that a' R a = .68659
#' ## vc = variance of composite
#' vc <- aRa
#' ## sf = scale factor
#' sf <- .68659/vc
#' a <- as.numeric(sqrt(sf)) * a
#' cat("\nKoopman Scaling\n")
#' print(round(a,2))
#' }
fungibleExtrema <- function(R.X, rxy, r.yhata.yhatb,  
                            Nstarts=100,  
                            MaxMin="Max", 
                            Seed = NULL,
                            maxGrad = 1E-05,
                            PrintLevel = 1){
 
  # Generate random seed if not supplied
  if(is.null(Seed)) Seed <- sample(1:10000, 1)
  set.seed(Seed)
  
  # Initialize convergenceStatus
  convergenceStatus = FALSE
  
  # auxiliary function definitions
  # vector norm
  norm <- function(x) x/as.numeric( sqrt(t(x) %*%x))
  # vector cosine 
  vec.cos <- function(x,y){ t(norm(x))%*%norm(y) }
  # vector length 
  lngth <- function(x) sqrt(t(x) %*%x)


##~~~~~~~~~~~~~Generate U matrix via Gram Schmidt~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
#    GenU.BAK <- function(mat,u){
#         ## 
#         p <- ncol(mat)
#         n <- nrow(mat)
#         oData <-matrix(0,n,p+1)
#         oData[,1]<-u
# 
#         for(i in 2:(p+1)){
#            oData[,i] <- resid(lm(mat[,(i-1)]~-1+oData[,1:(i-1)]))
#         }
# 
#         U<-oData[,2:(p+1)]
#         d <- diag(1/sqrt(diag(crossprod(U))))
#         U <- U%*%d
#         U
#     }#end GenU

##~~~~~~~~~~~~~Generate U matrix via QR Decomposition~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  GenU <- function(mat,u){
    qr.Q(qr(cbind(u,mat)))[,-1]
  }#end GenU

#~~~~~~~~~~~~~~~Compute Analytic Gradient~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

gradf <- function(sv){

   z<-sv[1:(NX-1)]       
   zpz <- t(z)%*%z 
   L <- sv[NX]        
   k <- r * u  + U %*% z * e     
   a <- V%*%Linv.sqrt%*%k
   apa <- t(a)%*%a
   ab <- as.numeric(t(a) %*%b.unit)
  
   dfdz <- (e * t(U) %*% Linv.sqrt %*%t(V)) %*% 
           (-as.numeric(ab*(apa)^-1.5) * a  + as.numeric((apa^-.5))*b.unit )
  
   
   if(MaxMin=="MIN") {
       dfdz <-   dfdz +4*L* as.numeric(zpz-1) * z  
       dfdL<-    (zpz - 1)^2  
   }       
   if(MaxMin=="MAX") {
    dfdz <-   dfdz -4 * L * as.numeric(zpz-1)*z  
    dfdL<-   -(zpz - 1)^2    
   } 
    
   c(dfdz,dfdL)      
 }   

##~~~~~~~~~~~~Main Function~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  #sv = start values
  minv <- function(sv){
      
      z<-sv[1:(NX-1)]
      
      zpz <- t(z)%*%z
      L <- sv[NX]     
      
      k <- r * u + + U %*% z * sqrt(1-r^2)
            
      k.star <- Linv.sqrt%*%k
      a <- V %*% k.star  
      len.a<-lngth(a)
      f <-  (as.numeric(len.a)^-1)*t(a)%*%b.unit  + L*(zpz - 1)^2 
      
      if(MaxMin=="MAX") f<- (as.numeric(len.a)^-1)*t(a)%*%b.unit  - L*(zpz - 1)^2 
      if(MaxMin=="MIN") f<- (as.numeric(len.a)^-1)*t(a)%*%b.unit  + L*(zpz - 1)^2  
      f
   
   }  ## end minv    
     
 

  MaxMin <- toupper(MaxMin)  
  NX <- ncol(R.X)

  #OLS weights
  b <- crossprod(solve(R.X),rxy)
  b.unit <- norm(b)
  
  len.b <- sqrt(t(b)%*%b)
  r <- as.numeric(r.yhata.yhatb)
  e <- sqrt(1-r^2)

  VLV <- eigen(R.X)
  V <- VLV$vectors
  L <- diag(VLV$values)
  Linv<-solve(L)
  Linv.sqrt <- solve( sqrt(L))
 
  u <- (sqrt(L)%*%t(V)%*%b)/ as.numeric(sqrt(t(b)%*%R.X%*%b))

  mat <- matrix(rnorm(NX*(NX-1)),NX,NX-1)
  U <- GenU(mat,u)
    
     FNSCALE <- 1
     if(MaxMin=="MAX")  FNSCALE <- FNSCALE * -1  
    
     # declare minf /maxf
     minf <- 99*FNSCALE
     maxf <- -1    
     breakflag<-0
     iter <- 0
  
    set.seed(Seed) 
    while(iter < Nstarts){
      # start values
      sv <- c(rnorm(NX-1)/sqrt(NX-1),1000)
      
      tmp <- try(optim(par=sv,
                       fn=minv,
                       gr=gradf,
                       method="BFGS",
                       control=list(fnscale=FNSCALE,
                                    maxit=500,
                                    parscale=c(rep(1,NX-1),1))))
                
                
                 if(abs(tmp$value)>1) tmp$convergence<-1
                 if((FNSCALE*tmp$value <= FNSCALE*minf) & (tmp$convergence==0))
                 { 
                    iter <- iter + 1
                    fdelta <-minf-tmp$value
                    minf<-tmp$value
                    out<-tmp           
                    z<-out$par[1:(NX-1)]
                    k <- r * u + U %*% z * sqrt(1-r^2)
                    a <- a.tilde <- V %*% Linv.sqrt %*% k 
                                       
                    scaling.weight <- (t(rxy) %*% a)/(t(a)%*%R.X%*%a)               
                    a <- as.numeric(scaling.weight) * a  
                    
                    max_grad = max(abs(gradf(out$par)))
                    
                    # PRINT
                    if(PrintLevel == 1){
                       cat(c(iter,"Current fnc val  = ", round(minf,5), "\n",
                       " Max abs gradient = ",  round(max_grad, 5),"\n"))
                    }
                    
                    # Check convergence status
                  
                    if( (max_grad <= maxGrad) && PrintLevel == 1){
                       convergenceStatus = TRUE
                       cat("\n ***Sucessful Convergence*** \n\n")
                       breakflag <- 1  
                    }  
                    
                 } # End tmp$minimum < minf
               
       if(breakflag) break
   }# End NLoop    
      

      if(sign(a[1])!=sign(a.tilde[1]))  a <- a*-1

      r.y.yhat <- sqrt( (t(b) %*% R.X %*%b) )
      

      
      
	  # compute vector cosine 
      s<-vcos(a,b) 
      ba.mat <- cbind(b,a)
      colnames(ba.mat) <- c("b","a")
      
    # compute final gradient 
      solution.gradient <- gradf(out$par)
      
      # PRINT
      if(PrintLevel == 1 && convergenceStatus == TRUE){
        
        cat("\n\n========= Results =============\n\n")
        if(MaxMin=="MAX") cat("  Maximizing cos(a,b):\n")
        if(MaxMin=="MIN") cat("  Minimizing cos(a,b):\n")        
        cat("  r(yhat.a,yhat.b) = ",round(r.yhata.yhatb,3),"\n")
        
        cat("  cos(a,b) = ",round(s,3),"\n")
        cat("  RSQb = ",round(r.y.yhat^2,3),"\n")
        cat("  RSQa = ",round((r.yhata.yhatb * r.y.yhat)^2,3),"\n")
        cat("  Relative loss = RSQb - RSQa = ",
             round( r.y.yhat^2-(r.yhata.yhatb * r.y.yhat)^2 ,3),"\n\n")
         
        print(ba.mat)
        
        cat("\n\n")
        cat("Analytic Gradient at Solution\n")
        print(matrix(solution.gradient,NX,1))
        cat("\n")
      } #END PRINT    
       
  # ---- RETURN ----
   list(cos.ab = s,
        a = a ,
        k = k,
        z = z,
        b = b,
        u = u,
        r.yhata.yhatb = r.yhata.yhatb,
        r.y.yhatb = r.y.yhat,
        gradient = solution.gradient,
        converged = convergenceStatus)
 } ## End of Fungible Function





# ##----------------------------------------##
# ##   EXAMPLE
# ##  This is Koopmnan's Table 2 Example
# ##----------------------------------------##
# 
# R.X <- matrix(c(1.00,  .69,  .49,  .39,
#                  .69, 1.00,  .38,  .19,
#                  .49,  .38, 1.00,  .27,
#                  .39,  .19,  .27, 1.00),4,4)
#                  
#                  
# b <- c(.39, .22, .02, .43)
# rxy <- R.X%*%b
# 
# OLSRSQ <- t(b)%*%R.X%*%b
# 
# #theta <- .02
# #r.yhata.yhatb <- sqrt( 1 - (theta)/OLSRSQ)
# 
# r.yhata.yhatb  <- .90
# output <- FungibleExtrema1(R.X, rxy, r.yhata.yhatb, Nstarts=500, MaxMin="Min")
# 
#  ## Scale to replicate Koopman
#       a<-output$a
#       a.old<-a
#       aRa<-t(a)%*%R.X%*%a
#       ##Scale a such that a' R a = .68659
#       ##vc =variance of composite
#       vc <- aRa
#       ##sf = scale factor
#       sf <- .68659/vc
#       a <- as.numeric(sqrt(sf))*a
#       cat("\nKoopman Scaling\n")
#       print(round(a,2))
# 
  
 


