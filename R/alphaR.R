# Niels Waller
# January 30, 2024
# 
#' Generate random R matrices with a known coefficient alpha
#'
#' @usage alphaR(alpha, k, Nmats, SEED)
#' 
#' @description
#' \code{alphaR} can generate a list of fungible correlation matrices with a 
#' user-defined  (standardized) coefficient \eqn{\alpha}. 
#' 
#' @param alpha  (numeric) A desired coefficient \eqn{\alpha} within 
#' the range \eqn{\alpha \in (-\infty, 1]}. 
#' @param k (integer). The order of each R (correlation) matrix. 
#' @param Nmats (integer) The number of fungible R matrices with a known 
#' \eqn{\alpha}. 
#' Default (\code{Nmats = 5}).
#' @param SEED (numeric)  The initial seed for the random number generator. 
#' If SEED is not supplied then the program will generate (and return) a randomly
#' generated seed.
#' 
#' @return
#'\itemize{
#'    \item  \strong{alpha} The desired (standardized) coefficient \eqn{\alpha}.
#'    \item \strong{R} The initial correlation matrix with a desired 
#'    coefficient \eqn{\alpha}. 
#'    \item \strong{Rlist} A list with \code{Nmats} fungible correlation 
#'    matrices with a desired coefficient \eqn{\alpha}.  
#'    \item  \strong{SEED} The initial value for the random number generator.
#' }  
#' @author Niels G. Waller
#'
#' @references  Waller, N. & Revelle, W. (2023). What are the mathematical 
#' bounds for coefficient \eqn{\alpha}? \emph{Psychological Methods}.
#'  doi.org/10.1037/met0000583
#' 
#' @examples
#'
#' ## Function to compute standardized alpha
#' Alphaz <- function(Rxx){
#'   k <- ncol(Rxx)
#'   k/(k-1) * (1 - (k/sum(Rxx)) ) 
#' }# END Alphaz
#'
#' ## Example 1
#' ## Generate 25 6 x 6 R matrices with a standardized alpha of .85
#' alpha =  .85   
#' k = 6
#' Nmats =  25 
#' SEED = 1
#'
#' out = alphaR(alpha, k , Nmats, SEED)
#' Alphaz(out$Rlist[[1]])
#'
#' ## Example 2
#' ## Generate 25 6 x 6 R matrices with a standardized alpha of -5
#' alpha =  -5   
#' k = 6
#' Nmats =  25 
#' SEED = 1
#'
#' out = alphaR(alpha, k , Nmats, SEED) 
#' Alphaz(out$Rlist[[5]])
#' 
#' @export
#' 
alphaR <- function(alpha = NULL, k = NULL, Nmats = 5, SEED = NULL){
  
  ## generate random seed if not supplied
  if(is.null(SEED)) SEED <- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) 
  set.seed(SEED)
  
  if(is.null(k)){
    stop("\n\nFATAL ERROR: k must be specified")
  }
  
  Rlist = NULL
  
# ---- Alphaz: Compute standardized alpha ----
  Alphaz <- function(Rxx){
    k <- ncol(Rxx)
    k/(k-1) * (1 - (k/sum(Rxx)) ) 
  }# END Alphaz

## Case I
if(alpha == 1){
  # R is a matrix of 1's
  Ralpha <- matrix(1, k, k) 
}

## Case II
if(alpha == 0 ){
  # R is an identity matrix
  Ralpha <- diag(k)
}

# Case III: alpha in (0,1) ----
if(alpha > 0 && alpha < 1){
  # S = sum(R) 
  S <- sumR <- k/(1 - alpha*(k-1)/k)
  r <- (S-k)/(k*(k-1))
  Ralpha <- matrix(r, k, k)
  diag(Ralpha) <- 1
}

# Case IV: alpha in (-infinity, 0) ----
if(alpha < 0){
  r <- (1 - k/(k-1))  # Common value of rij
# initial R
  Rx <- matrix(r, k, k)
  diag(Rx) <- 1

  c <- k/(2*(1 - alpha * (k-1)/k))
  # distribute 2c across all non diagonal r's
  eps <- 2*c/(k*(k-1))
  Ralpha <- Rx + eps; diag(Ralpha) <- 1
}

# Generate Nmats fungible R matrices with a known alpha
  if(Nmats > 1){
    Rlist = vector(mode="list", length = Nmats)
    for(i in 1:Nmats){
      Rlist[[i]] <-  fungible::Ravgr(Rseed = Ralpha, NVar = k)$R
    }
  }#END if(Nmats > 1)  
 
## RETURN VALUES
list( alpha =  Alphaz(Rxx = Ralpha),
      R = Ralpha,
      Rlist = Rlist,
      SEED = SEED)

} ## END alphaR

