##########################################################
# This function compute the asymptotic distribution-free #
# covariance matrix of covariances.                      #
#                                                        #
# Arguments:                                             #
# X - matrix of predictor scores                         #
# y - vector of criterion scores                         #
#                                                        #
# Output                                                 #
# adfCovMat - Asymptotic distribution-free estimate of   #
#             the covariance matrix                      #
##########################################################



#' Asymptotic Distribution-Free Covariance Matrix of Covariances
#' 
#' Function for computing an asymptotic distribution-free covariance matrix of
#' covariances.
#' 
#' 
#' @param X Data matrix.
#' @param y Optional vector of criterion scores.
#' @return \item{adfCovMat}{Asymptotic distribution-free estimate of the
#' covariance matrix of covariances}
#' @author Jeff Jones and Niels Waller
#' @references Browne, M. W. (1984). Asymptotically distribution-free methods
#' for the analysis of covariance structures. \emph{British Journal of
#' Mathematical and Statistical Psychology, 37,} 62--83.
#' @keywords Statistics
#' @examples
#' 
#' ## Generate non-normal data using monte1
#' set.seed(123)
#' 
#' ## we will simulate data for 1000 subjects
#' N <- 1000
#' 
#' ## R = the desired population correlation matrix among predictors
#' R <- matrix(c(1, .5, .5, 1), 2, 2)
#' 
#' ## Consider a regression model with coefficient of determination (Rsq):
#' Rsq <- .50
#' 
#' ## and vector of standardized regression coefficients
#' Beta <- sqrt(Rsq/t(sqrt(c(.5, .5))) %*% R %*% sqrt(c(.5, .5))) * sqrt(c(.5, .5))
#' 
#' ## generate non-normal data for the predictors (X)
#' ## x1 has expected skew = 1 and kurtosis = 3
#' ## x2 has expected skew = 2 and kurtosis = 5
#' X <- monte1(seed = 123, nvar = 2, nsub = N, cormat = R, skewvec = c(1, 2), 
#'            kurtvec = c(3, 5))$data
#'            
#' ## generate criterion scores 
#' y <- X %*% Beta + sqrt(1-Rsq)*rnorm(N)
#' 
#' ## Create ADF Covariance Matrix of Covariances
#' adfCov(X, y)
#' 
#' #>         11       12       13       22       23       33
#' #> 11 3.438760 2.317159 2.269080 2.442003 1.962584 1.688631
#' #> 12 2.317159 3.171722 2.278212 3.349173 2.692097 2.028701
#' #> 13 2.269080 2.278212 2.303659 2.395033 2.149316 2.106310
#' #> 22 2.442003 3.349173 2.395033 6.275088 4.086652 2.687647
#' #> 23 1.962584 2.692097 2.149316 4.086652 3.287088 2.501094
#' #> 33 1.688631 2.028701 2.106310 2.687647 2.501094 2.818664
#' @export
#' 
adfCov <- function(X, y=NULL) {
  Xy <- if(is.null(y)) X else cbind(X,y)
  dev <- scale(Xy,scale=FALSE)
  nvar <- ncol(dev)
  N <- nrow(dev)

# order of the covariance matrix of covariances
  ue <- nvar*(nvar + 1)/2

# container for indices
  s <- vector(length=ue, mode="character")

  z <- 0
  for(i in 1:nvar){
    for(j in i:nvar){
      z<-z+1
      s[z]<-paste(i,j,sep="")
    }
  }

# computes all possible combinations of the 
# indices in s
  v <- expand.grid(s, s)

# paste the index pairs togehter
  V <- paste(v[,1], v[,2], sep="")

# separate the indices into their own columns
  id.mat <- matrix(0,nrow=ue^2,4)
  for(i in 1:4) id.mat[,i] <- as.numeric(sapply(V,substr,i,i))

# create a matrix with the sequence of numbers 1:ue^2 by row
# we will use this to find the positions of the indices in id.mat
  M <- matrix(1:ue^2, ue, ue, byrow=TRUE)

# grabs the rows of the index pairs of interest
  r <- M[lower.tri(M, diag=TRUE)]

  ids <- id.mat[r,]

  adfCovMat <- adfCovBias <- matrix(0,ue,ue)
  covs <- covs.bias <- matrix(0,nrow(ids),1)

# compute the covariances using Browne (1984) Eqn 3.8

  for(i in 1:nrow(ids)) {
   
    w_ij <- crossprod(dev[,ids[i,1]],dev[,ids[i,2]])/N
    w_ik <- crossprod(dev[,ids[i,1]],dev[,ids[i,3]])/N
    w_il <- crossprod(dev[,ids[i,1]],dev[,ids[i,4]])/N
    w_jk <- crossprod(dev[,ids[i,2]],dev[,ids[i,3]])/N
    w_jl <- crossprod(dev[,ids[i,2]],dev[,ids[i,4]])/N
    w_kl <- crossprod(dev[,ids[i,3]],dev[,ids[i,4]])/N
     
  
    w_ijkl <- (t(dev[,ids[i,1]]*dev[,ids[i,2]])%*%
                 (dev[,ids[i,3]]*dev[,ids[i,4]])/N)         
    
    covs[i] <- (N*(N-1)*(1/((N-2)*(N-3)))*(w_ijkl - w_ij*w_kl) -
                  N*(1/((N-2)*(N-3)))*(w_ik*w_jl + w_il*w_jk - (2/(N-1))*w_ij*w_kl))
  }	

# Create Covariance Matrix

  adfCovMat[lower.tri(adfCovMat,diag=T)] <- covs	
  vars <- diag(adfCovMat)
  adfCovMat <- adfCovMat + t(adfCovMat) - diag(vars)

  ## add row and column labels
  rc<-expand.grid(1:nvar,1:nvar)
  rc<-unique(apply(t(apply(rc,1,sort)),1,paste, collapse=""))
  rc
  dimnames(adfCovMat)<-list(rc,rc)
  adfCovMat
  
}
