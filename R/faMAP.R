########################
#  faMAP
#  Niels Waller, June 16, 2009
# updated February 26, 2018
#         added recordPlot
#  Based on modified Matlab code orginally written by Brian O'Connor




#' Velicer's minimum partial correlation method for determining the number of
#' major components for a principal components analysis or a factor analysis
#' 
#' Uses Velicer's MAP (i.e., matrix of partial correlations) procedure to
#' determine the number of components from a matrix of partial correlations.
#' 
#' 
#' @param R input data in the form of a correlation matrix.
#' @param max.fac maximum number of dimensions to extract.
#' @param Print (logical) Print = TRUE will print complete results.
#' @param Plot (logical) Plot = TRUE will plot the MAP values.
#' @return \item{MAP}{Minimum partial correlations} \item{MAP4}{Minimum partial
#' correlations} \item{fm}{average of the squared partial correlations after
#' the first m components are partialed out.} \item{fm4}{see Velicer, Eaton, &
#' Fava, 2000.} \item{PlotAvgSq}{A saved object of the original MAP plot (based
#' on the average squared partial r's.)} \item{PlotAvg4th}{A saved object of
#' the revised MAP plot (based on the average 4th power of the partial r's.)}
#' @author Niels Waller
#' @references Velicer, W. (1976). Determining the number of components from
#' the matrix of partial correlations. Psychometrika, 41(3):321--327.
#' 
#' Velicer,W. F., Eaton, C. A. , & Fava, J. L. (2000). Construct explication
#' through factor or component analysis: A review and evaluation of alternative
#' procedures for determining the number of factors or components. In R. D.
#' Goffin & E. Helmes (Eds.). Problems and Solutions in Human Assessment:
#' Honoring Douglas N. Jackson at Seventy (pp. 41-71. Boston, MA: Kluwer
#' Academic.
#' @keywords statistics
#' @export
#' @importFrom graphics plot
#' @importFrom grDevices recordPlot
#' @examples
#' 
#' 	# Harman's data (1967, p 80) 
#' 	# R = matrix(c(
#' 	# 1.000,  .846,  .805,  .859,  .473,  .398,  .301,  .382,
#' 	#  .846, 1.000,  .881,  .826,  .376,  .326,  .277,  .415,
#' 	#  .805,  .881, 1.000,  .801,  .380,  .319,  .237,  .345,
#' 	#  .859,  .826,  .801, 1.000,  .436,  .329,  .327,  .365,
#' 	#  .473,  .376,  .380,  .436, 1.000,  .762,  .730,  .629,
#' 	#  .398,  .326,  .319,  .329,  .762, 1.000,  .583,  .577,
#' 	#  .301,  .277,  .237,  .327,  .730,  .583, 1.000,  .539,
#' 	#  .382,  .415,  .345,  .365,  .629,  .577,  .539, 1.000), 8,8)
#' 
#' 	  F <- matrix(c(  .4,  .1,  .0,
#' 	                  .5,  .0,  .1,
#' 	                  .6,  .03, .1,
#' 	                  .4, -.2,  .0,
#' 	                   0,  .6,  .1,
#' 	                  .1,  .7,  .2,
#' 	                  .3,  .7,  .1,
#' 	                   0,  .4,  .1,
#' 	                   0,   0,  .5,
#' 	                  .1, -.2,  .6, 
#' 	                  .1,  .2,  .7,
#' 	                 -.2,  .1,  .7),12,3)
#' 					 
#' 	  R <- F %*% t(F)
#' 	  diag(R) <- 1 
#' 	  
#'   	faMAP(R, max.fac = 8, Print = TRUE, Plot = TRUE) 
#' 
faMAP <- function(R, max.fac = 8, Print=TRUE, Plot=TRUE) {
	
nvars <- nrow(R)
ULU <- eigen(R)
eigval <- ULU$values
eigvect = ULU$vectors

I <- diag(nvars)

loadings = eigvect %*% diag(sqrt(eigval))
	
fm4 <- fm <- rep(0,max.fac)
fm[1]  <- sum(R^2 - I)/(nvars*(nvars-1))
fm4[1] <- sum(R^4 - I)/(nvars*(nvars-1))

for(m in 1:(max.fac-1)){
     biga <- loadings[,1:m]
     partcov = R - (biga %*% t(biga))
     d <- diag (  (1 / sqrt(diag(partcov))))
     pr <- d %*% partcov %*% d
     fm[m+1] <-  (sum(pr^2)-nvars)/(nvars*(nvars-1))
     fm4[m+1] <- (sum(pr^4)-nvars)/(nvars*(nvars-1))
 }	
 
 # identifying the smallest fm value & its location
 minfm.loc <- which.min(fm)
 minfm4.loc <- which.min(fm4)
 
 if(Print){
    cat("\nVelicer's Minimum Average Partial (MAP) Test\n\n")
    cat("The smallest average squared partial correlation is: ",
        round(min(fm),3),"\n")
    cat("The smallest average 4rth power partial correlation is: ",
        round(min(fm4),3),"\n\n")
 
    cat("The Number of Components According to the Original (1976) MAP Test is = ", minfm.loc,"\n")
    cat("The Number of Components According to the Revised  (2000) MAP Test is = ",minfm4.loc)
 }
 
 PlotAvgSq <- NULL
 
 m1 <- c("Original MAP (Avg squared partial r)", paste("\nNumber of Components = ", minfm.loc, sep=""))
 if(Plot){
 	plot(1:max.fac,fm,type="b", 
 	     main=m1,
 	     xlab="Dimensions",
 	     ylab="Avg squared partial r",
 	     xlim=c(1,max.fac))

  PlotAvgSq <- recordPlot()
 
 	m1 <- c("Revised MAP (Avg 4th partial r):\n", 
 	        paste("\nNumber of Components = ", 
 	        minfm4.loc, sep=""))
 		plot(1:max.fac,fm4,type="b", 
 	     main=m1,
 	     xlab="Dimensions",
 	     ylab="Avg 4th partial r") 
 		
 		PlotAvg4th <- recordPlot()
  }	
 
 invisible(list(MAP = minfm.loc, 
                MAP4 = minfm4.loc, 
                fm = fm, 
                fm4 = fm4,
                PlotAvgSq = PlotAvgSq,
                PlotAvg4th = PlotAvg4th))
}

