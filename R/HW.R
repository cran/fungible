#' Six data sets that  yield a Heywood case
#'
#' Six data sets that  yield a Heywood case in a 3-factor model. 
#'
#' @docType data
#' 
#' @usage data(HW)
#' 
#' @format Each data set is a  matrix with 150 rows and 12 variables:
#' \describe{ Each data set (HW1, HW2, ... HW6) represents a hypothetical sample
#' of 150 subjects from a population 3-factor model.  
#' The population factor loadings are given in \code{HW$popLoadings}.
#' }
#' 
#'
#' @keywords datasets
#
#' @examples 
#' data(HW)
#' 
#' # Compute a principal axis factor analysis 
#' # on the first data set  
#' RHW <- cor(HW$HW1)  
#' fapaOut <- faMain(R = RHW, 
#'                  numFactors = 3, 
#'                  facMethod = "fapa", 
#'                  rotate = "oblimin",
#'                  faControl = list(treatHeywood = FALSE))
#'
#'
#' fapaOut$faFit$Heywood
#' round(fapaOut$h2, 2)
#' 
#'
"HW"

