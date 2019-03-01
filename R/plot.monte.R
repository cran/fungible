#' Plot Method for Class Monte
#' 
#' plot method for class "monte"
#' 
#' 
#' @param x An object of class 'monte', usually, a result of a call to
#' \code{monte}.
#' @param \dots Optional arguments passed to plotting function.
#' @return The function \code{plot.monte} creates a scatter plot of matrices
#' plot (a splom plot).  Cluster membership is denoted by different colors in
#' the plot.
#' @import lattice
#' @export
#' @examples
#' 
#' #plot(monte.object)
#' 
plot.monte<-function(x,...){
#library(lattice)
z<-x   #monte.object
super.sym <- lattice::trellis.par.get("superpose.symbol")
grp.labs<-as.factor(paste("Group",z$data[,1],sep=" "))
lattice::splom(~z$data[,2:(z$nvar+1)], groups = grp.labs,  
      # panel = panel.superpose,
       key = list(title = "Monte: Scatter Plot Matrices",
                  columns = 3, 
                  points = list(pch = super.sym$pch[1:z$nclus],
                  col = super.sym$col[1:3]),text=list(as.vector(unlist(unique(grp.labs))))))

}
