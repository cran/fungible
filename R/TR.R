## TR
## Author: Niels Waller

#' @title Estimate the parameters of the Taylor-Russell function.    
#' 
#' @description
#' A Taylor-Russell function can be computed with any three of the following 
#' four variables: the Base Rate (BR); the Selection Ratio (SR); 
#' the Criterion Validity (CV) and the Positive Predictive Value (PPV). 
#' The \code{TR()} function will compute a Taylor Russell function 
#' when given any three of these parameters and estimate the remaining parameter. 
#'
#' @param BR (numeric): The Base Rate of successful criterion performance 
#'  (i.e., within the target population, the proportion of individuals who can 
#'  successfully execute the job demands). 
#' @param SR (numeric): The Selection Ratio. A real number between 0 and 1  
#'  that denotes the test selection ratio (i.e., the proportion of hired 
#'  candidates from the target population).
#' @param CV (numeric) The correlation (Criterion Validity) between the selection test 
#'  and a measure of job performance.
#' @param PPV (numeric): The Positive Predicted Value. The PPV denotes 
#'  the probability that a hired candidate has the necessary skills to 
#'  succeed on the job.
#' @param PrintLevel (integer): If \code{PrintLevel = 0} then no output will be printed to screen. 
#' If \code{PrintLevel = 1} then a brief summary of output is printed to screen. Default \code{PrintLevel = 1}. 
#' @param Digits (integer) Controls the number of significant digits 
#' in the printed output.   
#'  
#' @details  
#'  When any three of the main program arguments (BR, SR, CV, PPV) are specified (with the 
#'  remaining argument given a NULL value), \code{TR()} will calculate 
#'  the model-implied value for the remaining variable.  It will also compute the test Sensitivity 
#'  (defined as the probability that a qualified individual will be hired) and 
#'  test Specificity (defined as the probability that an unqualified individual 
#'  will not be hired), the True Positive rate, the False Positive rate, the 
#'  True Negative rate, and the False Negative rate.
#'   
#'  
#' @return
#'\itemize{
#'    \item  \strong{BR} The base rate.
#'    \item \strong{SR} The selection ratio.
#'    \item \strong{CV} The criterion validity.
#'    \item  \strong{PPV} The positive predictive value.
#'    \item \strong{Sensitivity} The test sensitivity rate. 
#'    \item \strong{Specificity}  The test specificity rate. 
#'    \item \strong{TP} The selection True Positive rate. 
#'    \item \strong{FP} The selection False Positive rate. 
#'    \item \strong{TN} The selection True Negative rate. 
#'    \item \strong{FN} The selection False Negative rate. 
#' }  
#' 
#'  
#' @author
#' \itemize{
#'   \item Niels G. Waller (nwaller@umn.edu)
#' }
#' 
#' @references
#' \itemize{
#'   \item Taylor, H. C. & Russell, J. (1939). The relationship of validity 
#'   coefficients to the practical effectiveness of tests in selection: 
#'   Discussion and tables. \emph{Journal of Applied Psychology, 23}, 565--578.
#'  }
#'  
#' @keywords stats 
#' @import mvtnorm
#' @export
#'  
#' @examples
#' 
#' ## Example 1:
#'        TR(BR = .3, 
#'           SR = NULL, 
#'           CV = .3, 
#'           PPV = .5,
#'           PrintLevel = 1,
#'           Digits = 3)
#'   
#' ## Example 2:
#' 
#'        TR(BR = NULL, 
#'           SR = .1012, 
#'           CV = .3, 
#'           PPV = .5,
#'           PrintLevel = 1,
#'           Digits = 3)
#'    
#' ## Example 3: A really bad test!
#'  # If the BR > PPV then the actual test
#'  # validity is zero. Thus, do not use the test!
#' 
#'        TR(BR = .50, 
#'           SR = NULL, 
#'           CV = .3, 
#'           PPV = .25,
#'           PrintLevel = 1,
#'           Digits = 3)   

TR <- function(BR = NULL,
                   SR = NULL, 
                   CV = NULL,
                   PPV = NULL,
                   PrintLevel = 1,
                   Digits = 3){
  
#  Enter NULL for any one variable, e.g.,
#   BR =  .2      # Base Rate       
#   SR  = .4      # Selection Ratio
#   CV   = NULL   # Criterion validity
#   PPV = .5      # Positive Predicted Value
#                  #  i.e., the proportion of positives 
#                  # that are True positives
  
 
#===============================================#  
# Taylor Russell functions
# ---- TR.sr ----  
TR.sr <- function(SR){
  threshold.x <- -qnorm(SR)
  threshold.y <- -qnorm(BR)


	TP <- mvtnorm::pmvnorm(lower = c(threshold.x,  threshold.y), 
	                      upper = c(+Inf, +Inf),
	                      sigma = matrix(c(1, CV, CV, 1),2,2))[1] 
	                      
	FP <- mvtnorm::pmvnorm(lower = c(threshold.x, -Inf), 
	                      upper = c(+Inf, threshold.y),
	                      sigma = matrix(c(1, CV, CV, 1),2,2))[1] 
	
	# function to minimize
	Q <- 500*( PPV - TP/(TP + FP) )^2  # return
  }# END TR.sr function

# ---- TR.br ----  
TR.br <- function(BR){
  
  threshold.x <- -qnorm(SR)
  threshold.y <- -qnorm(BR)
  
 
  TP <- mvtnorm::pmvnorm(lower = c(threshold.x,  threshold.y), 
                        upper = c(+Inf, +Inf),
                        sigma = matrix(c(1, CV, CV, 1),2,2))[1] 
  
  FP <- mvtnorm::pmvnorm(lower = c(threshold.x, -Inf), 
                        upper = c(+Inf, threshold.y),
                        sigma = matrix(c(1, CV, CV, 1),2,2))[1] 
  
  # function to minimize
  Q <- 500 * ( PPV - TP/(TP + FP) )^2  # return
  
}# END TR.br function

# ---- TR.r ----  
TR.r <- function(CV){
  
  threshold.x <- -qnorm(SR)
  threshold.y <- -qnorm(BR)
  
  
  TP <- mvtnorm::pmvnorm(lower = c(threshold.x,  threshold.y), 
                        upper = c(+Inf, +Inf),
                        sigma = matrix(c(1, CV, CV, 1),2,2))[1] 
  
  FP <- mvtnorm::pmvnorm(lower = c(threshold.x, -Inf), 
                        upper = c(+Inf, threshold.y),
                        sigma = matrix(c(1, CV, CV, 1),2,2))[1] 
  
  # function to minimize
  Q <- 500 * ( PPV - TP/(TP + FP) )^2  # return
  
}# END TR.r function


# ---- TR.ppv ----  
TR.ppv <- function(PPV){
  
  threshold.x <- -qnorm(SR)
  threshold.y <- -qnorm(BR)
  
  
  TP <- mvtnorm::pmvnorm(lower = c(threshold.x,  threshold.y), 
                        upper = c(+Inf, +Inf),
                        sigma = matrix(c(1, CV, CV, 1),2,2))[1] 
  
  FP <- mvtnorm::pmvnorm(lower = c(threshold.x, -Inf), 
                        upper = c(+Inf, threshold.y),
                        sigma = matrix(c(1, CV, CV, 1),2,2))[1] 
  
  # function to minimize
  Q <- 500 * ( PPV - TP/(TP + FP) )^2  # return
  
}# END TR.ppv function



# ---- Find Selection Ratio ----
 if(is.null(SR)){ 
     out <- optimize(f = TR.sr,
                     lower = 0.001,
                     upper = .999,
                     tol=1E-8,
                     maximum = FALSE)

     SR = out$minimum    
     
  } #END if SR == Null



# ---- Find Base Rate  ----
if(is.null(BR)){ 
  out <- optimize(f = TR.br, 
                  lower = 0.001, 
                  upper = .999,
                  maximum = FALSE)
  
  BR = out$minimum    
} #END if BR == Null

# ---- Find CV  ----
if(is.null(CV)){ 
  out <- optimize(f = TR.r, 
                  lower = -1, 
                  upper = 1,
                  maximum = FALSE)
  
  CV = out$minimum    
} #END if CV == Null

# ---- Find PPV  ----
if(is.null(PPV)){ 
  out <- optimize(f = TR.ppv, 
                  lower = 0, 
                  upper = 1,
                  maximum = FALSE)
  
  PPV = out$minimum    
} #END if PPV == Null


# ----Sensitivity and Specificity ----

 threshold.x <- -qnorm(SR)
 threshold.y <- -qnorm(BR)


  TP <- mvtnorm::pmvnorm(lower = c(threshold.x,  threshold.y), 
                       upper = c(+Inf, +Inf),
                       sigma = matrix(c(1, CV, CV, 1),2,2))[1] 

  FP <- mvtnorm::pmvnorm(lower = c(threshold.x, -Inf), 
                       upper = c(+Inf, threshold.y),
                       sigma = matrix(c(1, CV, CV, 1),2,2))[1] 

  Sensitivity = TP/BR
  Specificity <- (1-BR-FP) / (1-BR)
  
  FN = BR - TP
  TN = 1 - BR - FP

  
# ---- Error Checks ---- 
if(SR == 0 || BR == 0)
{ 
  
  stop("\n\n*** No solution in the parameter space. ***")
  
} 
if ( (BR > PPV) ){
  # The test has zero validity if the BR > PPV
  # These are the stats if you ignore the test
  # and select all applicants
  PPV = BR
  SR = 1
  CV = 0
  Sensitivity = 1
  Specificity = 0
  TP = NULL
  FP = NULL
  TN = NULL
  FN = NULL
  
  if(PrintLevel == 1){
    cat(
      "\n *** WARNING: BR > PPV ***\n",
      "Base Rate = ", round(BR,Digits), "\n",
      "Revised Positive Predicted Value  = ", round(PPV, Digits), "\n",
      "Selection Ratio = ", round(SR,Digits), "\n",  
      "Test validity = ZERO! Dont use the test!", "\n"
      )
  }#END if(PrintLevel == 1)
} #END if ( (BR > PPV) )

if ( (BR < PPV) & PrintLevel == 1 ){
  cat(
    "\n",
    "Base Rate = ", round(BR,Digits), "\n",
    "Selection Ratio = ", round(SR,Digits), "\n",  
    "Criterion Validity = ", round(CV,Digits), "\n",
    "Positive Predicted Value  = ", round(PPV, Digits),"\n",
    "Sensitivity = ", round(Sensitivity, Digits),"\n",
    "Specificity = ", round(Specificity, Digits),"\n",
    "TP = ", round(TP, Digits), "\n",
    "FP = ", round(FP, Digits), "\n",
    "TN = ", round(TN, Digits), "\n",
    "FN = ", round(FN, Digits))
    
}

invisible(list( 
      "BR" = BR,
      "SR" = SR,
      "CV" = CV,
      "PPV" = PPV,
      "Sensitivity" = Sensitivity,
      "Specifcity" = Specificity,
      "TP" = TP,
      "FP" = FP,
      "TN" = TN,
      "FN" = FN)
)

} ## END TR 
  

  
 
 