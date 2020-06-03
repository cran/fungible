#' Calculate Reference Eigenvalues for the Empirical Kaiser Criterion
#' 
#' 
#' @param R Input correlation matrix.
#' @param NSubj Number of subjects (observations) used to create R. 
#' @param Plot (logical). If \code{Plot = TRUE} the function will plot
#' the observed and reference eigenvalues of R.
#' @return 
#' \itemize{
#'  \item{ljEKC}, 
#'  \item{ljEKC1},
#'  \item{dimensions} The estimated number of common factors.
#'  }
#' @author Niels Waller
#' @keywords Statistics
#' @family Factor Analysis Routines
#' @references Braeken, J. & Van Assen, M. A.  (2017).  An empirical Kaiser criterion. 
#' \emph{Psychological Methods, 22}(3), 450-466.  
#' @export
#' @examples
#' 
#' data(AmzBoxes)
#' AmzBox20<- GenerateBoxData(XYZ = AmzBoxes[,2:4], 
#'                            BoxStudy = 20)$BoxData
#' RAmzBox20 <- cor(AmzBox20)
#' EKCout  <- faEKC(R = RAmzBox20, 
#'                 NSubj = 98,
#'                 Plot = TRUE)
#' 
#' 
faEKC <- function(R = NULL, NSubj = NULL, Plot = FALSE){
  
    # J = NVariables
    J <- nrow(R)
    N <- NSubj
    gamma <- J/N
  
    lup <- (1 + sqrt(gamma))^2

    l <- eigen(R)$values
    cumSumEig <- cumsum(l)
  
    ljRefFNC <-function(l, j){
      ((J - cumSumEig[j-1])/(J - j + 1)) * lup
    }  
    
    ljEKC1 <- ljEKC <- rep(99,J)
    for(i in 2:J){
      ljEKC1[i] <- max(ljRefFNC(l,i), 1)
      ljEKC[i] <- ljRefFNC(l,i)
    }  
    
    ljEKC[1] <- ljEKC1[1] <- lup
   
    dimensions <- sum(l > ljEKC1)
   
    # ----Plot ----
    if(Plot){
     plot(1:J, l,  
          type="b",
          ylab="Eigenvalues",
          xlab = "Dimensions",
          main = paste0("Empirical Kaiser Criterion\n ", dimensions, " Dimensions"))
      points(1:J, ljEKC1, type="b",col="red")
    }
   
   #----RETURN----
   list(ljEKC = ljEKC,
        ljEKC1 =  ljEKC1,
        dimensions = dimensions)
}## END faEKC





