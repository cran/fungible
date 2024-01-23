## Author: Niels Waller
## Version 1.1
## September 16, 2023
#' Generalized Taylor-Russell Function for Multiple Predictors
#' 
#'
#' @title A generalized (multiple predictor) Taylor-Russell function.  
#'
#' @param SR (vector)  A vector of Selection Ratios for N selection tests.
#' @param BR (scalar)  The Base Rate of criterion performance.
#' @param R  (matrix)  An (N + 1) x (N + 1) correlation matrix in which the 
#'  predictor/criterion correlations are in column  N + 1 of R.
#' @param PrintLevel (integer). If \code{PrintLevel = 0} then no output is 
#'  printed to screen. If \code{PrintLevel > 0} then output is printed to screen. 
#'  Defaults to \code{PrintLevel = 0}.
#' @param Digits (integer)  The number of significant digits in the printed 
#'  output. 
#'  
#' @return The following output variables are returned.
#'
#' \itemize{
#'   \item \strong{BR}: (scalar) The Base Rate of criterion performance. 
#'   \item \strong{SR}: (vector) The user-defined vector of predictor Selection 
#'   Ratios.
#'   \item \strong{R}: (matrix) The input correlation matrix.  
#'   \item \strong{TP}: (scalar) The percentage of True Positives.
#'   \item \strong{FP}: (scalar) The percentage of False Positives.
#'   \item \strong{TN}: (scalar) The percentage of True Negatives.
#'   \item \strong{FN}: (scalar) The percentage of False Negatives.
#'   \item \strong{Accepted}: The percentage of selected individuals 
#'   (i.e., TP + FP).
#'   \item \strong{PPV}: The Positive Predictive Value. This is the probability 
#'   that a selected individual is a True Positive.
#'   \item \strong{Sensitivity}: The test battery Sensitivity rate.  This is the 
#'   probability that a person who is acceptable on the criterion is called  
#'   acceptable by the test battery.
#'   \item \strong{Specificity}: The test battery Specificity rate.  This is the 
#'   probability that a person who falls below the criterion threshold 
#'   is deemed unacceptable by the test battery.
#' }
#' 
#' @author
#' \itemize{
#'   \item Niels G. Waller (nwaller@umn.edu)
#' }
#' @references
#' \itemize{
#'   \item Taylor, H. C. & Russell, J. (1939). The relationship of validity 
#'   coefficients to the practical effectiveness of tests in selection: 
#'   Discussion and tables. Journal of Applied Psychology, 23(5), 565--578.
#'   
#'   \item Thomas, J. G., Owen, D., & Gunst, R. (1977). Improving the use 
#'   of educational tests as selection tools. Journal of Educational 
#'   Statistics, 2(1), 55--77. 
#' }
#' @keywords stats
#' @import mvtnorm
#' @export
#' 
#' @examples
#' # Example 1
#' # Reproduce Table 3 (p. 574) of Taylor and Russell
#' 
#' r <- seq(0, 1, by = .05)
#' sr <- c(.05, seq(.10, .90, by = .10), .95)
#' num.r <- length(r)
#' num.sr <- length(sr)
#' 
#' old <- options(width = 132)
#' 
#' Table3 <- matrix(0, num.r, num.sr)
#' for(i in 1 : num.r){
#'    for(j in 1:num.sr){
#'    
#'      Table3[i,j] <-  TaylorRussell(
#'                        SR = sr[j],
#'                        BR = .20, 
#'                        R = matrix(c(1, r[i], r[i], 1), 2, 2), 
#'                        PrintLevel = 0,
#'                        Digits = 3)$PPV  
#'    
#'   }# END over j
#' }# END over i
#'
#' rownames(Table3) <- r
#' colnames(Table3) <- sr
#' Table3 |> round(2)
#'
#' # Example 2
#' # Thomas, Owen, & Gunst (1977) -- Example 1: Criterion = GPA
#' 
#' R <- matrix(c(1, .5, .7,
#'              .5, 1, .7,
#'             .7, .7, 1), 3, 3)
#'
#'  # See Table 6: Target Acceptance = 20%
#'  out.20 <- TaylorRussell(
#'  SR = c(.354, .354),  # the marginal probabilities
#'  BR = .60, 
#'  R = R,
#'  PrintLevel = 1) 
#'
#' # See Table 6:  Target Acceptance = 50%
#' out.50 <- TaylorRussell(
#'  SR = c(.653, .653),   # the marginal probabilities
#'  BR = .60, 
#'  R = R,
#'  PrintLevel = 1) 
#'  
#'  options(old)
#'  
 TaylorRussell <- function(
                 SR = NULL,
                 BR = NULL, 
                 R = NULL, 
                 PrintLevel = 0,
                 Digits = 3){  
  
  # to get consistent results from repeated function calls,
  # initialize the seed
  # see the pmvnorm {mvtnorm} help page 
  # set.seed(66)
   
  # Convert selection ratios to thresholds
  threshold.x <- stats::qnorm(SR, lower.tail = FALSE)
  threshold.y <- stats::qnorm(BR, lower.tail = FALSE)
  
  Num.x <- length(threshold.x)
#===============================================#  
# Taylor Russell functions

  
  TP <- mvtnorm::pmvnorm(lower = 
                          c(threshold.x,  threshold.y), 
                          upper =  rep(+Inf, Num.x + 1),
                          mean = rep(0, Num.x +1 ),
                          corr = R,
                          keepAttr=FALSE)
  
  FP <- mvtnorm::pmvnorm(lower = c(threshold.x, 
                                  -Inf), 
                        upper = c(rep(+Inf, Num.x), 
                                  threshold.y),
                        mean = rep(0, Num.x +1 ),
                        corr  = R,
                        keepAttr=FALSE)
  
# ---- PPV ---- 
# Probability of being a TP if called a Positive  
  PPV = TP/(TP + FP)


# ---- Sensitivity ---- 
  # what proportion of positives are called positive
  Sen <- TP/BR

# ---- Specificity ---- 
  # what proportion of negatives are called negative
  # Proportion of Negatives = 1-BR
  Spe <- (1 - BR - FP)/(1 - BR)

# ---- Accepted ---- 
  Accepted <- TP + FP
  
  FN = BR - TP
  TN = 1 - TP - FP - FN
  
if(PrintLevel > 0){
  BR.prnt <- round(BR, Digits)
  SR.prnt <- round(SR, Digits)
  PPV.prnt <- round(PPV, Digits)
  Sen.prnt <- round(Sen, Digits)
  Spe.prnt <- round(Spe, Digits)
  

     cat(
       " BR   =        ", BR.prnt, "\n",
       "SR    =       ",  SR.prnt, "\n",
       "TP    =       ", round(TP,Digits),"\n",
       "FP    =       ", round(FP,Digits),"\n",
       "TN    =       ", round(TN,Digits),"\n",
       "FN    =       ", round(FN,Digits),"\n",
       "Accepted =    ", round(Accepted, Digits),"\n",
       "PPV  =        ", PPV.prnt, "\n",
       "Sensitivity = ", Sen.prnt, "\n",
       "Specificity = ", Spe.prnt, "\n")
}  
  
invisible( list(
    "BR" = BR,
    "SR" = SR,  
    "R" = R,
    "TP" = TP,
    "FP" = FP,
    "TN" = TN,
    "FN" = FN,
    "Accepted" = Accepted,
    "PPV" = PPV,
    "Sensitivity" = Sen,
    "Specificity" = Spe))

}# END TaylorRussell function


