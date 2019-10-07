#' Factor Pattern Matrix for Thurstone's 26  box attributes. 
#'
#' Factor Pattern Matrix for Thurstone's 26  box attributes. 
#'
#' @docType data
#' 
#' @usage data(ThurstoneBox26)
#' 
#' @format The original factor pattern (3 graphically rotated centroid factors) from Thurstone's 26 hypothetical box data as 
#' reported by Thurstone (1947, p. 371). The so-called Thurstone invariant box problem 
#' contains measurements on the following 26 functions of length (x), width (y), and height (z). 
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
#' 
#' @keywords datasets
#' 
#' @details 
#' Two data sets have been described in the literature as Thurstone's Box Data 
#' (or Thurstone's Box Problem). The first consists of 20 measurements on a set of 20 
#' hypothetical boxes (i.e., Thurstone made up the data).  Those data are available 
#' in \strong{Box20}. The second data set was collected by 
#' Thurstone to provide an illustration of the invariance of simple structure 
#' factor loadings. In his classic textbook on multiple factor analysis 
#' (Thurstone, 1947), Thurstone states that ``[m]easurements of a random collection 
#' of thirty boxes were actually made in the Psychometric Laboratory and recorded 
#' for this numerical example. The three dimensions, x, y, and z, were recorded 
#' for each box. A list of 26 arbitrary score functions was then prepared'' (p. 369). The 
#' raw data for this example were not published.  Rather, Thurstone reported a 
#' correlation matrix for the 26 score functions (Thurstone, 1947, p. 370). Note that, presumably 
#' due to rounding error in the reported correlations, the correlation matrix 
#' for this example is non positive definite. This file includes the rotated centroid solution 
#' that is reported in his book (Thurstone, 1947, p. 371).
#' 
#' @references 
#' Thurstone, L. L.  (1947).  Multiple factor analysis.  Chicago: University of Chicago Press. 
#' 
#' @seealso \code{\link{Box20}}, \code{\link{AmzBoxes}}
#' 
#' @examples 
#' data(ThurstoneBox26)  
#' ThurstoneBox26
"ThurstoneBox26"


