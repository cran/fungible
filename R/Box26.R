#' R matrix for Thurstone's 26 hypothetical box attributes. 
#'
#' Correlation matrix for Thurstone's 26 hypothetical box attributes. 
#'
#' @docType data
#' 
#' @usage data(Box26)
#' 
#' @format Correlation matrix for Thurstone's 26 hypothetical box attributes. 
#' The so-called Thurstone invariant box problem contains measurements on the 
#' following 26 functions of length, width, and height. 
#'   \strong{Box26} variables:
#'   \enumerate{
#'   \item x
#'   \item y
#'   \item z
#'   \item xy 
#'   \item xz
#'   \item yz 
#'   \item x^2 * y
#'   \item x * y^2
#'   \item x^2 * z
#'   \item x * z^ 2
#'   \item y^2 * z
#'   \item y * z^2
#'   \item x/y
#'   \item y/x
#'   \item x/z
#'   \item  z/x
#'   \item  y/z
#'   \item  z/y
#'   \item 2x + 2y
#'   \item 2x + 2z
#'   \item 2y + 2z
#'   \item sqrt(x^2 + y^2)
#'   \item sqrt(x^2 + z^2)
#'   \item sqrt(y^2 + z^2)
#'   \item xyz
#'   \item sqrt(x^2 + y^2 + z^2)
#'   }
#' \itemize{
#'    \item \strong{x} Box length
#'    \item \strong{y} Box width
#'    \item \strong{z} Box height
#'  }  
#' 
#'
#' @keywords datasets
#
#' @examples 
#' 
#' data(Box26)
#' fout <- faMain(R     = Box26,
#'                numFactors    = 3,
#'                facMethod     = "fapa",
#'                rotate        = "varimax",
#'                bootstrapSE   = FALSE,
#'         rotateControl = list(
#'                numberStarts = 100,  
#'                standardize  = "none"),
#'                Seed = 123)
#'
#'         print( round(fout$loadings, 2 ) ) 
#'         
#' fout <- faMain(R     = Box26,
#'                numFactors    = 3,
#'                facMethod     = "fapa",
#'                rotate        = "varimax",
#'                bootstrapSE   = FALSE,
#'         rotateControl = list(
#'                numberStarts = 100,  
#'                standardize  = "CM"),
#'                Seed = 123)
#'
#'         print( round(fout$loadings, 2 ) ) 
#' 
#'
"Box26"

