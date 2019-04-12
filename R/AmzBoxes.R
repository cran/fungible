#' Length, width, and height measurements for 98 Amazon shipping boxes
#'
#' Length, width, and height measurements for 98 Amazon shipping boxes
#'
#' @docType data
#' 
#' @usage data(AmzBoxes)
#' 
#' @format A data set of measurements for 98 Amazon shipping boxes. 
#' These data were downloaded from the  BoxDimensions website: (\url{https://www.boxdimensions.com/}).
#' The data set includes five variables:
#' \itemize{
#'    \item Amazon Box Size
#'    \item Length (inches) 
#'    \item Width (inches) 
#'    \item Height (inches)
#'    \item Volume (inches)
#'  }  
#' 
#'
#' @keywords datasets
#
#' @examples 
#' data(AmzBoxes)
#' 
#' hist(AmzBoxes$`Length (inches)`,
#'      main = "Histogram of Box Lengths",
#'      xlab = "Length",
#'      col = "blue")
#' 
#'
"AmzBoxes"

