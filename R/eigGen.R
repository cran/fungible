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


