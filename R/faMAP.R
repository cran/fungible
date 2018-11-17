
########################
#  faMAP
#  Niels Waller, June 16, 2009
# updated February 26, 2018
#         added recordPlot
#  Based on modified Matlab code orginally written by Brian O'Connor


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

