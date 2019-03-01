#---------------------------------------------------------------#
# restScore  Author Niels G Waller
# Latest revision January 27, 2016
#---------------------------------------------------------------#


#' Plot an ERF using rest scores
#' 
#' Plot an empirical response function using rest scores.
#' 
#' 
#' @param data N(subjects)-by-p(items) matrix of 0/1 item response data.
#' @param item Generate a rest score plot for item \code{item}.
#' @param NCuts Divide the rest scores into \code{NCuts} bins of equal width.
#' @return A restscore plot with 95\% confidence interval bars for the
#' conditional probability estimates. \item{item}{The item number.}
#' \item{bins}{A vector of bin limits and bin sample sizes.} \item{binProb}{A
#' vector of bin conditional probabilities.}
#' @author Niels Waller
#' @keywords statistics
#' @export
#' @importFrom graphics abline plot axis lines points
#' @importFrom methods is
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
#' data<-genFMPData(NSubj = NSubj, bParam = b, seed = 345)$data
#' 
#' ## generate a rest score plot for item 12.
#' ## the grey horizontal lines in the plot
#' ## respresent pseudo asymptotes that
#' ## are significantly different from the 
#' ## (0,1) boundaries
#' restScore(data, item = 12, NCuts = 9)
#' 
restScore<- function(data,  item, NCuts=10){

	Nitems<-ncol(data)
   tot.score<-apply(data,1,sum)
   pmat<-rep(0,NCuts)
   rest.score<-tot.score - data[,item]
 
    probs<-unlist(lapply(split(data[,item],cut(rest.score,NCuts)),mean))
    group.N<-unlist(lapply(split(data[,item],cut(rest.score,NCuts)),length))

    se.p<- sqrt(( probs * (1-probs))/group.N) 
    plot(1:NCuts,probs,axes=FALSE,ylim=c(0,1),
        xlab="Rest Scores",
        ylab="Probability  of  a  Keyed  Response",
        pch=16,col="blue",cex=1.5,
        lwd=3,
        main=paste("Item ",item, sep="")) 

    # bound limits of std error bars to (0,1)
    upperProbs <- probs + 1.96*se.p
    upperProbs[upperProbs>1]<-1
    lowerProbs <- probs - 1.96*se.p
    lowerProbs[lowerProbs < 0] <- 0

    for (i in 1:NCuts) { 
        lines( c((1:NCuts)[i],(1:NCuts)[i]),c(probs[i],upperProbs[i]),lwd=3, lty=4) 
        }
    for (i in 1:NCuts) { 
        lines( c((1:NCuts)[i],(1:NCuts)[i]),c(probs[i],lowerProbs[i]),lwd=3, lty=4) 
        }         

    
    axis(side=1,at=1:NCuts, pos=0,
      labels=names(table(cut(rest.score,NCuts))),
      cex.axis=.68,
      lwd=3)
    
    
    axis(side=2,lwd=3,pos= .75,cex.axis=.68)
    abline(h=1,lwd=3)
    abline(h=0,lwd=3)

    points(1:NCuts,probs,type="l",lwd=3, pch=1, 
      cex=1.5,
      ylim=c(0,1))
   
   # draw lines at upper and lower probability bounds
    if((probs[1]-1.96*se.p[1])>0)
     abline(h=probs[1]-1.96*se.p[1],col="lightgrey",lty=3,lwd=3)
    if((probs[NCuts]+1.96*se.p[NCuts])<1)
     abline(h=probs[NCuts]+1.96*se.p[NCuts],col="grey",lty=3,lwd=3)
  
    invisible(list(item = item, bins = table(cut(rest.score,NCuts)),binProb = probs ))
}    


