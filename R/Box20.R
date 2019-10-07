#' Length, width, and height measurements for Thurstone's 20 boxes
#'
#' Length, width, and height measurements for Thurstone's 20 hypothetical boxes
#'
#' @docType data
#' 
#' @usage data(Box20)
#' 
#' @format A data set of measurements for Thurstone's 20 hypothetical boxes. 
#' The data set includes three variables:
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
#' data(Box20)
#' 
#' hist(Box20$x,
#'      main = "Histogram of Box Lengths",
#'      xlab = "Length",
#'      col = "blue")
#' 
#' # To create the raw data for Thurstone's 20 hypothetical 
#' # box attributes:
#' data(Box20)
#'  ThurstoneBox20 <- GenerateBoxData(XYZ = Box20,
#'                                  BoxStudy = 20,
#'                                  Reliability = 1,
#'                                  ModApproxErrVar = 0)$BoxData  
#'
#' RThurstoneBox20 <- cor(ThurstoneBox20)   
#'
#' # Smooth matrix to calculate factor indeterminacy values
#' RsmThurstoneBox20 <- smoothBY(RThurstoneBox20)$RBY
#' 
#' fout <- faMain(R = RsmThurstoneBox20,
#'               numFactors = 3,
#'               rotate = "varimax",
#'               facMethod = "faregLS",
#'               rotateControl = list(numberStarts = 100,
#'                                    maxItr =15000))
#' summary(fout, digits=3)
#' 
#' # Note that given the small ratio of subjects to variables,
#' # it is not possible to generate data for this example with model error 
#' # (unless SampleSize is increased).
"Box20"

