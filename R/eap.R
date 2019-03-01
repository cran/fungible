#' Compute eap trait estimates for FMP and FUP models
#' 
#' Compute eap trait estimates for items fit by filtered monotonic polynomial
#' IRT models.
#' 
#' 
#' @param data N(subjects)-by-p(items) matrix of 0/1 item response data.
#' @param bParams A p-by-9 matrix of FMP or FUP item parameters and model
#' designations. Columns 1 - 8 hold the (possibly zero valued) polynomial
#' coefficients; column 9 holds the value of \code{k}.
#' @param NQuad Number of quadrature points used to calculate the eap
#' estimates.
#' @param priorVar Variance of the normal prior for the eap estimates. The
#' prior mean equals 0.
#' @param mintheta,maxtheta NQuad quadrature points will be evenly spaced
#' between \code{mintheta} and \code{maxtheta}
#' @return \item{eap trait estimates.}{}
#' @author Niels Waller
#' @keywords statistics
#' @export
#' @examples
#' 
#' 
#' ## this example demonstrates how to calculate 
#' ## eap trait estimates for a scale composed of items 
#' ## that have been fit to FMP models of different 
#' ## degree 
#' 
#' NSubjects <- 2000
#' 
#' ## Assume that 
#' ## items 1 - 5 fit a k=0 model,
#' ## items 6 - 10 fit a k=1 model, and 
#' ## items 11 - 15 fit a k=2 model.
#' 
#' 
#'  itmParameters <- matrix(c(
#'   #  b0    b1     b2    b3    b4  b5, b6, b7,  k
#'   -1.05, 1.63,  0.00, 0.00, 0.00,  0,     0,  0,   0, #1
#'   -1.97, 1.75,  0.00, 0.00, 0.00,  0,     0,  0,   0, #2
#'   -1.77, 1.82,  0.00, 0.00, 0.00,  0,     0,  0,   0, #3
#'   -4.76, 2.67,  0.00, 0.00, 0.00,  0,     0,  0,   0, #4
#'   -2.15, 1.93,  0.00, 0.00, 0.00,  0,     0,  0,   0, #5
#'   -1.25, 1.17, -0.25, 0.12, 0.00,  0,     0,  0,   1, #6
#'    1.65, 0.01,  0.02, 0.03, 0.00,  0,     0,  0,   1, #7
#'   -2.99, 1.64,  0.17, 0.03, 0.00,  0,     0,  0,   1, #8
#'   -3.22, 2.40, -0.12, 0.10, 0.00,  0,     0,  0,   1, #9
#'   -0.75, 1.09, -0.39, 0.31, 0.00,  0,     0,  0,   1, #10
#'   -1.21, 9.07,  1.20,-0.01,-0.01,  0.01,  0,  0,   2, #11
#'   -1.92, 1.55, -0.17, 0.50,-0.01,  0.01,  0,  0,   2, #12
#'   -1.76, 1.29, -0.13, 1.60,-0.01,  0.01,  0,  0,   2, #13
#'   -2.32, 1.40,  0.55, 0.05,-0.01,  0.01,  0,  0,   2, #14
#'   -1.24, 2.48, -0.65, 0.60,-0.01,  0.01,  0,  0,   2),#15
#'   15, 9, byrow=TRUE)
#'  
#' # generate data using the above item parameters
#' ex1.data<-genFMPData(NSubj = NSubjects, bParams = itmParameters, 
#'                     seed = 345)$data
#' 
#' 
#' ## calculate eap estimates for mixed models
#' thetaEAP<-eap(data = ex1.data, bParams = itmParameters, 
#'                    NQuad = 25, priorVar = 2, 
#'                    mintheta = -4, maxtheta = 4)
#' 
#' ## compare eap estimates with initial theta surrogates
#' 
#' if(FALSE){     #set to TRUE to see plot
#' 
#'   thetaInit <- svdNorm(ex1.data)
#'   plot(thetaInit,thetaEAP, xlim = c(-3.5,3.5), 
#'                          ylim = c(-3.5,3.5),
#'                          xlab = "Initial theta surrogates",
#'                          ylab = "EAP trait estimates (Mixed models)")
#' }                         
#' 
eap<-function(data, bParams, NQuad = 21, priorVar =2, mintheta=-4, maxtheta=4){
  
  ##  Vers January 25, 2016
  ## eapMixed can be used to calculate eap theta estimates for item sets that fit 
  ##          monotonic polynomial IRT models of possibly different degree
  ##
  ## bParams an NItems by 9 matrix. Cols 1 : 8 include  coefficients of
  ##        the polynomial
  ##
  ## bParams col 9 should contain k (polynomial of degree 2k+1)

  if(is.data.frame(bParams)) bParams <- as.matrix(bParams)
  
  if(ncol(bParams)!=9) stop("\n\nbParams should have 9 columhns")
  
  breakPoints<-seq(mintheta,maxtheta,by=(maxtheta-mintheta)/NQuad)
  centers<-breakPoints+(maxtheta-mintheta)/(2*NQuad)
  #remove last center to create:
  #Quadrature points
  QuadPt<-centers[-length(centers)]
  
  
  NItem <- ncol(data)
  NSubj <- nrow(data)
  
  # inter QuadPt distance
  Quadinterval<-QuadPt[2] - QuadPt[1]
  # quadrature weights
  QuadWght <- dnorm(QuadPt, mean=0, sd=sqrt(priorVar))
  QuadWght <- QuadWght/sum(QuadWght)
  
  theta <- QuadPt
  thetaMat <- matrix(c(rep(1,NQuad), theta, theta^2, theta^3, theta^4, theta^5, theta^6, theta^7), nrow=NQuad, ncol=8)
  
  
  P <- function(m){
    1/(1+exp(-m))
  }
  
  eap <- vector(mode="numeric",length = NSubj)
  
  #initialize at 1
  EAPVec0 <- rep(1, length=NQuad)
  
  for(iSubj in 1:NSubj){
 
    EAPVec <- EAPVec0
   
    for(item in 1:NItem){
       mVec <- thetaMat %*% bParams[item,1:8]
       PVec <-P(mVec)   
       EAPVec <- EAPVec * ( PVec^data[iSubj,item]  *  (1-PVec)^(1 - data[iSubj,item ]))
    }
    eap[iSubj]<-sum(EAPVec *QuadWght * QuadPt)/sum(EAPVec * QuadWght)
  }
  
  #return eap estimates
  eap
  
}

