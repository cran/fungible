###################################################################
#---------------------------------------------------------------#
# empirical Response Function
#---------------------------------------------------------------#


#' Utility fnc to compute the components for an empirical response function
#' 
#' Utility function to compute empirical response functions.
#' 
#' 
#' @param theta Vector of estimated latent trait scores.
#' @param data A matrix of binary item responses.
#' @param whichItem Data for an erf will be generated for whichItem.
#' @param min Default = -3. Minimum value of theta.
#' @param max Default = 3. Maximum value of theta.
#' @param Ncuts Number of score groups for erf.
#' @return \item{probs}{A vector (of length Ncuts) of bin response
#' probabilities for the empirical response function.} \item{centers}{A vector
#' of bin centers. } \item{Ni}{Bin sample sizes.} \item{se.p}{Standard errors
#' of the estimated bin response probabilities.}
#' @author Niels Waller
#' @keywords statistics
#' @export
#' @examples
#' 
#' NSubj <- 2000
#' 
#' #generate sample k=1 FMP  data
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
#' theta <- rnorm(NSubj)  
#' data<-genFMPData(NSubj = NSubj, bParam = b, theta = theta, seed = 345)$data
#' 
#' erfItem1 <- erf(theta, data, whichItem = 1, min = -3, max = 3, Ncuts = 12)
#' 
#' plot( erfItem1$centers, erfItem1$probs, type="b", 
#'       main="Empirical Response Function",
#'       xlab = expression(theta),
#'       ylab="Probability",
#'       cex.lab=1.5)
#' 
erf<- function(theta,  data, whichItem, min=-3, max=3, Ncuts=12){
  
  breakPoints<-seq(min,max,by=(max-min)/Ncuts)
  centers<-breakPoints+(max-min)/(2*Ncuts)
  # remove last center
  centers<-centers[-length(centers)]
  
  probs<-unlist(lapply(split(data[,whichItem],cut(theta,breaks=breakPoints)),mean))
  group.N<-unlist(lapply(split(data[,whichItem],cut(theta,breaks=breakPoints)),length))
  se.p<- sqrt(( probs * (1-probs))/group.N) 
  #     for (i in 1:Ncuts) { 
  #         lines( c((1:Ncuts)[i],(1:Ncuts)[i]),c(probs[i],(probs[i]+1.96*se.p[i]) ),lwd=line.width) 
  #         }
  #     for (i in 1:Ncuts) { 
  #         lines( c((1:Ncuts)[i],(1:Ncuts)[i]),c(probs[i],(probs[i]-1.96*se.p[i]) ),lwd=line.width) 
  #         }         
  list(probs=probs, centers = centers, Ni=group.N, se.p = se.p)
}    


