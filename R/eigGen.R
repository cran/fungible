##########################################################################
## eigGen:  generate eigenvalues for R matrices with underlying component structure
## Niels Waller
## updated November 28, 2017 (eig of 1st Minor factor must be < 1)
## October 9, 2016
##
## nDimensions     : total number of dimensions (variables)
## nMajorFactors   : number of major factors 
## PrcntMajor      : percent of variance accounted for by major factors
## threshold       : minimm difference in eigenvalues between last major factor
##                 : and first minor factor


#' Generate eigenvalues for R matrices with underlying component structure
#' 
#' Generate eigenvalues for R matrices with underlying component structure
#' 
#' 
#' @param nDimensions Total number of dimensions (variables).
#' @param nMajorFactors Number of major factors.
#' @param PrcntMajor Percentage of variance accounted for by major factors.
#' @param threshold Minimm difference in eigenvalues between the last major
#' factor and the first minor factor.
#' @return A vector of eigenvalues that satisfies the above criteria.
#' @author Niels Waller
#' @keywords Statistics
#' @export
#' @examples
#' 
#' ## Example
#' set.seed(323)
#' nDim <- 25   # number of dimensions
#' nMaj <- 5    # number of major components
#' pmaj <- 0.70 # percentage of variance accounted for
#'              # by major components
#' thresh <- 1  # eigenvalue difference between last major component 
#'              # and first minor component
#'  
#' L <- eigGen(nDimensions = nDim, nMajorFactors = nMaj, 
#'             PrcntMajor = pmaj, threshold = thresh)
#' 
#' maxy <- max(L+1)
#' 
#' plotTitle <- paste("  n Dimensions = ", nDim, 
#'                    ",  n Major Factors = ", nMaj, 
#' 				           "\n % Variance Major Factors = ", pmaj*100, 
#' 						   "%", sep = "")
#' 				 
#' plot(1:length(L), L, 
#'      type = "b", 
#'      main = plotTitle,
#'      ylim = c(0, maxy),
#'      xlab = "Dimensions", 
#' 	   ylab = "Eigenvalues",
#' 	   cex.main = .9)				 
#' 
#' 
eigGen<-function(nDimensions = 15, nMajorFactors = 5, 
                   PrcntMajor =.80, threshold =.5){
  
  nMinorFactors <- nDimensions - nMajorFactors                
  delta <- 0
  eigs <- rep(99, nDimensions)
  
  while( (delta < threshold)  |  (eigs[nMajorFactors + 1] >= 1) ) {
    
    eigMajor <-sort(runif(nMajorFactors,0,1),decreasing = TRUE)
    eigMajor <-eigMajor/sum(eigMajor)
    
    eigMinor <- sort(runif(nMinorFactors,0,1),decreasing = TRUE)
    eigMinor <- eigMinor/sum(eigMinor)
    
    eigs <- sort(nDimensions*c(PrcntMajor*eigMajor, (1-PrcntMajor)*eigMinor),decreasing = TRUE)
    
    delta <- eigs[nMajorFactors] - eigs[nMajorFactors+1]
   
   
  }

  eigs
}  
##########################################################################


