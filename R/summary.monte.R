#' Summary Method for an Object of Class Monte
#' 
#' summary method for class ``monte"
#' 
#' 
#' @aliases summary.monte print.summary.monte
#' @param object An object of class \code{monte}, usually, a result of a call
#' to \code{monte}.
#' @param digits Number of digits to print. Default = 3.
#' @param compute.validities Logical: If TRUE then the program will calculate
#' the indicator validities (eta^2) for the generated data.
#' @param Total.stats Logical: If TRUE then the program will return the
#' following statistics for the total sample: (1) indicator correlation matrix,
#' (2) indicator skewness, (3) indicator kurtosis.
#' @param \dots Optional arguments.
#' @return Various descriptive statistics will be computed within groups
#' including"
#'  
#' @return 
#'  \describe{
#'   \item{clus.size}{Number of objects within each group.}
#'   \item{centroids}{Group centroids.}
#'   \item{var.matrix}{Within group variances.}
#'   \item{cor.list}{Expected within group correlations.}
#'   \item{obs.cor}{Observed within group correlations.}
#'   \item{skew.list}{Expected within group indicator skewness values.}
#'   \item{obs.skew}{Observed within group indicator skewness values.}
#'   \item{kurt.list}{Expected within group indicator kurtosis values.}
#'   \item{obs.kurt}{Observed within group indicator kurtosis values.}
#'   \item{validities}{Observed indicator validities.}
#'   \item{Total.cor}{Total sample correlation matrix.}
#'   \item{Total.skew}{Total sample indicator skewness.}
#'   \item{Total.kurt}{Total sample indicator kurtosis.}
#'  }
#' @export
#' 
#' @examples
#' 
#' ## set up a 'monte' run for the Fisher iris data
#' 
#' sk.lst <- list(c(0.120,  0.041,  0.106,  1.254),                     #
#'                 c(0.105, -0.363, -0.607, -0.031),
#'                 c(0.118,  0.366,  0.549, -0.129) )
#'               
#'               
#' kt.lst <- list(c(-0.253, 0.955,  1.022,  1.719),
#'                 c(-0.533,-0.366,  0.048, -0.410),
#'                 c( 0.033, 0.706, -0.154, -0.602))
#'                 
#' cormat <- lapply(split(iris[,1:4],iris[,5]), cor)
#' 
#' my.iris <- monte(seed = 123, nvar = 4, nclus = 3, cor.list = cormat, 
#' 	              clus.size = c(50, 50, 50),
#'                 eta2 = c(0.619, 0.401, 0.941, 0.929), 
#'                 random.cor = FALSE,
#'                 skew.list = sk.lst, kurt.list = kt.lst, 
#'                 secor = .3, 
#'                 compactness = c(1, 1, 1), 
#'                 sortMeans = TRUE)
#' summary(my.iris)
#' 
summary.monte<-function(object,digits=3,compute.validities=FALSE,Total.stats=TRUE, ...){
z<-object    #monte.object
nclus<-z$nclus
nvar<-z$nvar

clus<-z$data[,1]
obs.cor.list<-as.list(rep(0,nclus))
obs.skew.matrix<-obs.kurt.matrix<-matrix(0,nclus,nvar)                                  
                

for(i in 1:nclus){  
  obs.cor.list[[i]]<-cor(z$data[clus==i,2:(nvar+1)])
  obs.skew.matrix[i,]<-apply(z$data[clus==i,2:(nvar+1)],2,skew)
  obs.kurt.matrix[i,]<-apply(z$data[clus==i,2:(nvar+1)],2,kurt)
}

ans <- z["call"]
ans$obs.cor<-obs.cor.list
ans$obs.skew<-obs.skew.matrix
ans$obs.kurt<-obs.kurt.matrix

cat("\nCall monte:","\n")
print(z$call)

cat("\nSeed = ",z$seed,"\n")

cat("\n\nNumber of objects in each group:")
for(i in 1:nclus){
 cat("\nGroup", i,"=",z$clus.size[i])
 }


##---------------------------------------------
##          GROUP CENTROIDS
##---------------------------------------------

centroids <-matrix(99,nclus,nvar)
cat("\n\n\nGroup Centroids\n")
for(i in 1:nclus){
  cat("Group ",i, "  ")
  centroids[i,] <- apply(z$data[ z$data[,1]==i, -1],2,mean)
  cat(round(centroids[i,],digits),"\n")
}


##---------------------------------------------
##          WITHIN GROUP VARIANCES
##---------------------------------------------

var.matrix<-matrix(99,nclus,nvar)
cat("\n\nWithin Group Variances\n")
for(i in 1:nclus){
  cat("Group ",i, "  ")
  var.matrix[i,]<- round(apply(z$data[ z$data[,1]==i, -1],2,var),digits)
  cat(var.matrix[i,],"\n")
}

RowNames <- paste0("Group ", 1:nclus)
ColNames <- paste0("var", 1:nvar)
dimnames(var.matrix) <- list(RowNames, ColNames)

grp.pairs <- expand.grid(1:z$nclus,1:z$nclus)
grp.pairs <- grp.pairs[grp.pairs[,1]!=grp.pairs[,2],]
grp.pairs <- grp.pairs[sort.list(grp.pairs[,1]),]

cat("\n\nRatio of Within Group Variances")
for(i in 1:length(grp.pairs[,1])){
   cat("\n\nGroup ",grp.pairs[i,1], "/ Group", grp.pairs[i,2], "  ")
   cat(round( var.matrix[ grp.pairs[i,1], ]/ var.matrix[ grp.pairs[i,2], ], digits))
}   



cat("\n\n\nExpected within group correlations:","\n")

for(i in 1:nclus){
  cat("\nGroup ",i,"\n")
  if(is.null(z$cor.list)) 
    cat("Identity Matrix\n")
  else
  print(round(z$cor.list[[i]],digits))
}


cat("\n\nObserved within group correlations:","\n")
for(i in 1:nclus){
  cat("\nGroup ",i)
  print(round(obs.cor.list[[i]],digits))
}

cat("\nExpected within group indicator skewness:","\n")
  if(is.null(z$skew.list))
    cat("All expected skew = 0\n")
  else {
    ex.sk.mat<-matrix(unlist(z$skew.list),nclus,nvar,byrow=TRUE)    
    dimnames(ex.sk.mat)<-list(paste("Group",1:nclus,sep=" "),paste("V",1:nvar,sep=""))
    print(ex.sk.mat)
  }  


cat("\nObserved within group indicator skewness:","\n")
dimnames(obs.skew.matrix)<-list(paste("Group",1:nclus,sep=" "),paste("V",1:nvar,sep=""))
dimnames(obs.kurt.matrix)<-list(paste("Group",1:nclus,sep=" "),paste("V",1:nvar,sep=""))
print(round(obs.skew.matrix,digits))

cat("\n\nExpected within group indicator kurtosis:","\n")
  if(is.null(z$kurt.list))
     cat("All expected kurtosis = 0\n")
  else{ 
    ex.kt.mat<-matrix(unlist(z$kurt.list),nclus,nvar,byrow=TRUE)
    dimnames(ex.kt.mat)<-list(paste("Group",1:nclus,sep=" "),paste("V",1:nvar,sep=""))
    print(ex.kt.mat)
  }  

cat("\nObserved within group indicator kurtosis:","\n")
print(round(obs.kurt.matrix,digits))



##----------------------------------------------------------
##                   VALIDITIES
##----------------------------------------------------------

ans$validities<- "Not computed"
if(compute.validities){
   validities <- rep(0,z$nvar)
     for(i in 1:z$nvar){
         validities[i]<-summary(lm(z$data[,i+1]~as.factor(z$data[,1])))$r.sq
     }

     cat("\nObserved indicator validities:","\n")
     print(round(validities,digits))
     ans$validities<-validities
     ans$eta2 <- z$eta2
 }  
 
 
 
 if(Total.stats){
 ##--------------------------------------------------------
 ##            Total Sample Results
 ##------------------------------------------------------
 cat("\n\nStats for Total Sample\n")
 cat("\nTotal sample correlation matrix")
 ans$Total.cor<-cor( z$data[,-1] )
 print(ans$Total.cor,digits)
 
 cat("\nTotal sample indicator skewness\n")
 ans$Total.skew<-apply(z$data[,-1],2,skew)
 cat( round(ans$Total.skew,digits) )
 
 
 cat("\n\nTotal sample indicator kurtosis\n")
  ans$Total.kurt<-apply(z$data[,-1],2,kurt)
 cat(round(ans$Total.kurt,digits) )
 cat("\n\n")
 }
 
 ans$seed<-z$seed
 ans$clus.size<-z$clus.size
 ans$centroids<-centroids
 ans$var.matrix <- var.matrix
 ans$cor.list<-z$cor.list
 ans$expected.skew.list <-z$ex.sk.mat
 ans$expected.kurtosis.list <-z$ex.kt.mat
 ans$skew.list<-z$skew.list
 ans$kurt.list<-z$kurt.list
 ans$Total.stats<-Total.stats

 
 
 
 class(ans)<-"summary.monte"
invisible(ans)
}
