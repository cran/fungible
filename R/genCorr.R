#################################################
# Generating correlation matrices               #
# with prespecified Eigen Values                #
# Based on Marsaglia and Olkin (1984) Algorithm #
#                                               #
# Jeff Jones                                    #
#                                               #
# Input:                                        #
#                                               #
#   eigenval - vector of desired eigen values   #
#   seed     - set.seed(seed)                   #
#                                               #
# Output:                                       #
#                                               #
#	 R        - correlation matrix               #
#################################################
#################################################




#' Generate Correlation Matrices with User-Defined Eigenvalues
#' 
#' Uses the Marsaglia and Olkin (1984) algorithm to generate correlation
#' matrices with user-defined eigenvalues.
#' 
#' 
#' @param eigenval A vector of eigenvalues that must sum to the order of the
#' desired correlation matrix.  For example: if you want a correlation matrix
#' of order 4, then you need 4 eigenvalues that sum to 4. A warning message
#' will display if sum(eigenval) != length(eigenval)
#' @param seed Either a user supplied seed for the random number generator or
#' `rand' for a function generated seed. Default seed=`rand'.
#' @return Returns a correlation matrix with the eigen-stucture specified by
#' eigenval.
#' @author Jeff Jones
#' @references Jones, J. A. (2010). GenCorr: An R routine to generate
#' correlation matrices from a user-defined eigenvalue structure. \emph{Applied
#' Psychological Measurement, 34}, 68-69.
#' 
#' Marsaglia, G., & Olkin, I. (1984). Generating correlation matrices.
#' \emph{SIAM J. Sci. and Stat. Comput., 5}, 470-475.
#' @keywords datagen
#' @export
#' @examples
#' 
#' 
#' ## Example
#' ## Generate a correlation matrix with user-specified eigenvalues
#' set.seed(123)
#' R <- genCorr(c(2.5, 1, 1, .3, .2))
#' 
#' print(round(R, 2))
#' 
#' #>       [,1]  [,2]  [,3]  [,4]  [,5]
#' #> [1,]  1.00  0.08 -0.07 -0.07  0.00
#' #> [2,]  0.08  1.00  0.00 -0.60  0.53
#' #> [3,] -0.07  0.00  1.00  0.51 -0.45
#' #> [4,] -0.07 -0.60  0.51  1.00 -0.75
#' #> [5,]  0.00  0.53 -0.45 -0.75  1.00
#' 
#' print(eigen(R)$values)
#' 
#' #[1] 2.5 1.0 1.0 0.3 0.2
#' 
#' 
genCorr <- function(eigenval,seed="rand"){
	
  if(!isTRUE(all.equal(sum(eigenval),length(eigenval)))) stop("Sum of eigenvalues not equal to Number of Variables\n")
  if(seed!="rand") set.seed(seed)

  norm <- function(x) x/as.numeric( sqrt(x %*%t(x)))

## step (i) in Marsaglia & Olkin

  Nvar <- length(eigenval)
  L <- diag(eigenval)
  E <- diag(Nvar)
  P <- matrix(0,Nvar,Nvar)

  if(identical(eigenval,rep(1,Nvar))) return(diag(Nvar))
  else {
  for(i in 1:(Nvar-1)) {

      xi <- as.vector(rnorm(Nvar))                  # step (ii) 
      xi <- xi%*%E                                  # create vector in rowspace of E

      a <- (xi %*% (diag(Nvar)-L) %*% t(xi))

      d_sq <- -1

      while(d_sq <= 0) {                             # if b^2 - ac < 0, recompute eta

        eta <- as.vector(rnorm(Nvar))                # step (iii)
        eta <- eta%*%E                               # create another vector in rowspace of E

        b <- (xi %*% (diag(Nvar)-L) %*% t(eta))      # step (iv)

        cc <- (eta %*% (diag(Nvar)-L) %*% t(eta))

        d_sq <- ((b^2) - (a*cc))

      } # end while

    ## step (v)

    r <- as.numeric(((b + sign(runif(1,-1,1))*sqrt(d_sq)))/(a))
    zeta.tmp <- (r*xi - eta)
    zeta <- sign(runif(1,-1,1))*norm(zeta.tmp)

    P[i,] <- zeta
    E <- E - t(zeta)%*%zeta
  
  } # end for

  P[Nvar,] <- norm(as.vector(rnorm(Nvar))%*%E)      # normalized random vector in rowspace of E

  R <- P%*%L%*%t(P)                                 # correlation matrix

  R     
  }

}

