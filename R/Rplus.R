# Niels Waller
# March 25, 2024

#' Generate a random R matrix with all rij > 0
#' 
#' Rplus(Seed, NVar = NULL,  rdist =  NULL,  param1 = NULL,  
#'       param2 = NULL, s = NULL)
#'       
#'
#' @param Seed (scalar) A seed for the random number generator.  If 
#' \code{Seed = NULL} then the program will generate a random seed. 
#' @param NVar (integer)  The order of the generated correlation (\strong{R}) matrix. 
#' @param rdist (character) A one or two letter character string that determines
#' the distribution for the elements of the Cholesky factor of the
#' randomly generated covariance matrix. Possible values include
#' \itemize{
#'  \item \strong{B} Beta(param1, param2)
#'  \item \strong{X} Chisq(param1
#'  \item \strong{E} Exponential(param1)
#'  \item \strong{F} F(param1, param2)
#'  \item \strong{G} Gamma(param1, param2)
#'  \item \strong{IB} Inverse Beta(param1, param2)
#'  \item \strong{IX} Inverse Chisq(param1)
#'  \item \strong{LN} Log Normal(param1, param2)
#'  \item \strong{U} Uniform(param1, param2)
#'  \item  \strong{W} Weibull(param1, param2)
#' }
#' @param param1 (scalar). The first parameter of \strong{rdist}.
#' @param param2 (scalar). The second parameter of \strong{rdist}
#' @param  s (scalar). An optional scalar \eqn{s \in [0,1]} that
#' scales the generated correlations by a factor of \strong{s}. 
#' Default (\code{s = NULL}). If \code{s} is supplied then no 
#' correlations in \strong{R} will be larger than \code{s}.
#' 
#' @return 
#' \itemize{
#'  \item \strong{R}  The generated \strong{R} matrix.  
#'  \item \strong{NVar} The order of \strong{R}.
#'  \item \strong{rdist} The user-supplied character strong for \strong{rdist},
#'  \item \strong{param1} The first parameter of \strong{rdist}. 
#'  \item \strong{param2} The possible second parameter of \strong{rdist}.
#'  \item \strong{s} The user-supplied scaling factor.
#'  \item \strong{Seed} the starting value for the random number generator.
#' }
#' 
#' @author Niels G. Waller
#'
#' @examples
#'  # --- Example 1: Generate two random R with positive rij ----
#'  R = Rplus(Seed = 123,  
#'            NVar = 6,  
#'            rdist =  "F",  
#'            param1 = 2.5,  
#'            param2 = 25,
#'            s = .8)$R
#'
#'  R |> round(2)
#' 
#' R = Rplus(Seed = 456,  
#'           NVar = 6,  
#'           rdist =  "F",  
#'           param1 = 2.5,  
#'           param2 = 25,
#'           s = .8)$R
#'        
#'  R |> round(2)   
#'
#' @importFrom LaplacesDemon rinvbeta rinvgamma rinvchisq            
#' @export           



#------------------------
Rplus <- function(Seed, NVar = NULL,  rdist =  NULL,  
                  param1 = NULL,  param2 = NULL,
                  s = NULL){
  
    if(is.null(Seed)) Seed <- sample(1:1E6, 1)
    set.seed(Seed)
    
    if( !(rdist %in% c("B", "X", "E", "F",
                     "G", "IB", "IG", "IX",
                     "LN", "U", "W"))  ){
      stop("\n\n*** Invalid entry for 'rdist'\n")
    }
  
    
    K <- matrix(0, NVar, NVar)
    # Number of independent elements of R
    d = NVar * (NVar + 1)/2
    
    
    # ---- rdist choices -------
    # Create a randomly generate Cholesky factor (K) of a covariance matrix.
    # The lower triangle of K is populated with realizations from
    # one of the following distributions.
    
    # ---- Beta ----    
    if(rdist == "B" && !is.null(param1) && !is.null(param2)){
      K[lower.tri(K, diag=TRUE)] <- rbeta(d, param1, param2)
    }
    
    # ---- Chisq ----
    if( rdist == "X" && !is.null(param1)  ){
      K[lower.tri(K, diag=TRUE)] <- rchisq(d, param1)
    } 
    
    
    # ---- Exponential ----    
    if(rdist == "E" && !is.null(param1) ){
      K[lower.tri(K, diag=TRUE)] <- rexp(d, param1)
    }
    
    # ---- F ---
    if(rdist == "F" && !is.null(param1) && !is.null(param2)){
      x <- rf(d, param1, param2)
      K[lower.tri(K, diag=TRUE)] <- x
    }
    
    
    # ---- Gamma ----    
    if(rdist == "G" && !is.null(param1) && !is.null(param2)){
      K[lower.tri(K, diag=TRUE)] <- rgamma(d, param1, param2)
    }
    
    
    # ---- Inverse Beta ---
    if(rdist == "IB" && !is.null(param1) && !is.null(param2)){
      K[lower.tri(K, diag=TRUE)] <- LaplacesDemon::rinvbeta(d, param1, param2)
    }
    
    # ---- Inverse Gamma ---
    if(rdist == "IG" && !is.null(param1) && !is.null(param2)){
      K[lower.tri(K, diag=TRUE)] <- LaplacesDemon::rinvgamma(d, param1, param2)
    }
    
    # ---- Inverse Chisq ---
    if(rdist == "IX" && !is.null(param1)){
      K[lower.tri(K, diag=TRUE)] <- LaplacesDemon::rinvchisq(d, param1)
    }
    
    # ---- Log Normal ----    
    if(rdist == "LN" && !is.null(param1) && !is.null(param2)){
      K[lower.tri(K, diag=TRUE)] <-  rlnorm(d, param1, param2)
    }
    
    # ---- Uniform ----
    if( rdist == "U" && !is.null(param1) && !is.null(param2) ){
      K[lower.tri(K, diag=TRUE)] <- runif(d, param1, param2)
    }
    
    # ---- Weibull ----    
    if(rdist == "W" && !is.null(param1) && !is.null(param2)){
      K[lower.tri(K, diag=TRUE)] <- rweibull(d, param1, param2)
    }
    
    
    # Generate a covariance matrix from K and then convert it into a cor matrix
    R <- cov2cor(K %*% t(K))
    
    # Permute rows and cols of R so that all rij have the
    # same marginal distributions
    P <- diag(NVar)
    r <- sample(1:NVar, NVar, replace = FALSE)
    P <- P[r,] #Permutation matrix
    R <- P %*% R %*% t(P)
    
    # Shrink rij so that max(rij) <= s
    if(!is.null(s)){
      R <- s * R
      diag(R) <- 1
    }
    
    
    # Return
    invisible(list(R = R,
                   NVar = NVar,
                   rdist = rdist,
                   param1 = param1,
                   param2 = param2,
                   s = s,
                   Seed = Seed))
} ## END  Rplus   


