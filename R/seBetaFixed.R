#------------------------------------------------------#
# This function computes the Cov of std Reg coeffs for #
# fixed X  using a method described by Yuan and Chan   #
# (2011)                                               #
#                                                      #
#                                                      #
# Arguments:                                           # 
# X            - matrix of predictor scores            #  
# y            - vector of criterion scores            #   
# cov.x        - covariance matrix for predictors      #    
# cov.xy       - vector of covariances between         #  
#                predictors and criterion              #  
# var.y        - criterion variance                    #  
# var.error    - var.error                             #
# Nobs         - number of observations                # 
#                                                      #
# This function accepts either (1) raw data, or (2)    #
# second-order moments (covariances or correlations)   #
# and sample size.                                     #
#                                                      #
# Output                                               #
#  seBeta      - normal theory standard errors for     #
#                standarized regression coefficients   #
#                with Fixed predictors.                #
#  covBeta     - normal theory covariance matrix of    #
#                standarized regsion coefficients      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





#' Covariance Matrix and Standard Errors for Standardized Regression
#' Coefficients for Fixed Predictors
#' 
#' Computes Normal Theory Covariance Matrix and Standard Errors for
#' Standardized Regression Coefficients for Fixed Predictors
#' 
#' 
#' @param X Matrix of predictor scores.
#' @param y Vector of criterion scores.
#' @param cov.x Covariance or correlation matrix of predictors.
#' @param cov.xy Vector of covariances or correlations between predictors and
#' criterion.
#' @param var.y Criterion variance.
#' @param var.error Optional argument to supply the error variance: var(y -
#' yhat).
#' @param Nobs Number of observations.
#' @return \item{cov.Beta}{Normal theory covariance matrix of standardized
#' regression coefficients for fixed predictors.} \item{se.Beta}{Standard
#' errors for standardized regression coefficients for fixed predictors.}
#' @author Jeff Jones and Niels Waller
#' @seealso \code{\link{seBeta}}
#' @references Yuan, K. & Chan, W. (2011). Biases and standard errors of
#' standardized regression coefficients. \emph{Psychometrika, 76(4)}, 670-690.
#' @keywords Statistics
#' @export
#' @import MASS
#' @examples
#' 
#' ## We will generate some data and pretend that the Predictors are being held fixed
#' 
#' library(MASS)
#' R <- matrix(.5, 3, 3); diag(R) <- 1
#' Beta <- c(.2, .3, .4)
#' 
#' rm(list = ".Random.seed", envir = globalenv()); set.seed(123)
#' X <- mvrnorm(n = 200, mu = rep(0, 3), Sigma = R, empirical = TRUE)
#' y <- X %*% Beta + .64*scale(rnorm(200))
#' 
#' seBetaFixed(X, y)
#' 
#' # $covBeta
#' #              b1           b2           b3
#' # b1  0.003275127 -0.001235665 -0.001274303
#' # b2 -0.001235665  0.003037100 -0.001491736
#' # b3 -0.001274303 -0.001491736  0.002830157
#' # 
#' # $seBeta
#' #         b1         b2         b3 
#' # 0.05722872 0.05510989 0.05319922
#' 
#' ## you can also supply covariances instead of raw data
#' 
#' seBetaFixed(cov.x = cov(X), cov.xy = cov(X, y), var.y = var(y), Nobs = 200)
#' 
#' # $covBeta
#' #              b1           b2           b3
#' # b1  0.003275127 -0.001235665 -0.001274303
#' # b2 -0.001235665  0.003037100 -0.001491736
#' # b3 -0.001274303 -0.001491736  0.002830157
#' # 
#' # $seBeta
#' #         b1         b2         b3 
#' # 0.05722872 0.05510989 0.05319922
#' 
#' 
seBetaFixed <- function(X=NULL,y=NULL,cov.x=NULL,cov.xy=NULL,
	                    var.y=NULL,var.error = NULL,Nobs=NULL) {

########################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~ Error Checking ~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  if(is.null(X) & !is.null(y)) 
    stop("\n y is not defined\n Need to specify both X and y\n")
  
  if(!is.null(X) & is.null(y)) 
    stop("\n X is not defined\n Need to specify both X and y\n")
			
  if(is.null(X) & 	is.null(y)) {
	
    if(is.null(cov.x) | is.null(cov.xy) | is.null(var.y) | is.null(Nobs))
      stop("\nYou need to specify all covariances and the number of observations\n")
      
    N <- Nobs	
    p <- nrow(cov.x)
		
  } else {
    
# if function is given data, compute all the covariances and variances    
	
    cov.x <- cov(X)
    cov.xy <- cov(X,y)
    var.y <- var(y)
    N <- length(y)
    p <- ncol(X)

  }

# create diagonal matrix of predictor standard deviations (Dx),
# unstandardized regression coefficients (b), and
# error variance from regression model (var.error)

  var.y <- as.numeric(var.y)
  Dx <- diag(sqrt(diag(cov.x)))	

  b <- solve(cov.x)%*%cov.xy
  
# if error variance is not supplied, calculate it
  
  if(is.null(var.error)) var.error <- as.numeric(var.y - t(b)%*%cov.x%*%b)*(N-1)/(N-p-1)	

# create Gamma matrix which is the 
# covariance matrix of b (unstandardized) and sigma^2 y (variance of y)

  var.b <- solve(cov.x)*var.error
  cov.b.var.y <- 2*b*var.error  
  var.var.y <- 4*(t(b)%*%cov.x%*%b)*var.error + 2*var.error^2
  Gamma <- rbind(cbind(var.b,cov.b.var.y),c(cov.b.var.y,var.var.y))

# Create jacobian for the function g(b,var.y) = Dx%*%b/sd.y

  deriv_b <- Dx*(1/sqrt(var.y))
  deriv_var.y <- -Dx%*%b/(2*var.y^(3/2))	
	
  jacobian <- cbind(deriv_b, deriv_var.y)	
	
# create covariance matrix of standardized regression coefficients
  
  cov.mat <- jacobian%*%Gamma%*%t(jacobian)/N

# name the covariance matrix  
  
  b.nms <- NULL
  for(i in 1:p) b.nms[i] <- paste("b",i,sep='')

  rownames(cov.mat) <- colnames(cov.mat) <- b.nms
  
  list(covBeta=cov.mat, seBeta=sqrt(diag(cov.mat)))	
	
}	

