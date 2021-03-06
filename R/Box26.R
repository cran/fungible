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
#' @details 
#' Two data sets have been described in the literature as Thurstone's Box Data 
#' (or Thurstone's Box Problem). The first consists of 20 measurements on a set of 20 
#' hypothetical boxes (i.e., Thurstone made up the data).  Those data are available 
#' in \strong{Box20}. The second data set,  which is described in this help file, was collected by 
#' Thurstone to provide an illustration of the invariance of simple structure 
#' factor loadings. In his classic textbook on multiple factor analysis 
#' (Thurstone, 1947), Thurstone states that ``[m]easurements of a random collection 
#' of thirty boxes were actually made in the Psychometric Laboratory and recorded 
#' for this numerical example. The three dimensions, x, y, and z, were recorded 
#' for each box. A list of 26 arbitrary score functions was then prepared'' (p. 369). The 
#' raw data for this example were not published.  Rather, Thurstone reported a 
#' correlation matrix for the 26 score functions (Thurstone, 1947, p. 370). Note that, presumably 
#' due to rounding error in the reported correlations, the correlation matrix 
#' for this example is non positive definite.
#' 
#' @references 
#' Thurstone, L. L.  (1947).  Multiple factor analysis.  Chicago: University of Chicago Press. 
#' @keywords datasets
#' 
#' @seealso \code{\link{Box20}}, \code{\link{AmzBoxes}}
#' @family Factor Analysis Routines
#' 
#' @examples 
#' 
#' data(Box26)
#' fout <- faMain(R     = Box26,
#'                numFactors    = 3,
#'                facMethod     = "faregLS",
#'                rotate        = "varimax",
#'                bootstrapSE   = FALSE,
#'         rotateControl = list(
#'                numberStarts = 100,  
#'                standardize  = "none"),
#'                Seed = 123)
#'
#' summary(fout)  
#'     
#' # We now choose Cureton-Mulaik row standardization to reveal 
#' # the underlying factor structure. 
#'           
#' fout <- faMain(R     = Box26,
#'                numFactors    = 3,
#'                facMethod     = "faregLS",
#'                rotate        = "varimax",
#'                bootstrapSE   = FALSE,
#'         rotateControl = list(
#'                numberStarts = 100,  
#'                standardize  = "CM"),
#'                Seed = 123)
#'
#' summary(fout)  
#' 
#'
"Box26"

